from __future__ import print_function
from srfpython.standalone.database import Database
from srfpython.depthdisp.depthmodels import depthmodel_from_arrays
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
import warnings
import numpy as np
import time, os


# -------------------------------------
# run file : based on sqlite
# -------------------------------------
class RunFile(Database):

    def drop(self):
        """erase existing tables and recreate them"""

        self.cursor.execute('drop table if exists PARAMVALUES')
        self.cursor.execute('drop table if exists DISPVALUES')
        self.cursor.execute('drop table if exists PARAMETERS')
        self.cursor.execute('drop table if exists DISPPOINTS')
        self.cursor.execute('drop table if exists MODELS')
        self.cursor.execute('drop table if exists CHAINS')

        self.cursor.execute('''create table DISPPOINTS (
            POINTID      integer primary key autoincrement,
            WAVE         varchar not null, --R/L
            TYPE         varchar not null, --C/U
            MODE         int not null,  
            FREQ         real not null, --Hz
            UNIT         varchar not null, --velocity
            constraint U unique (WAVE, TYPE, MODE, FREQ),
            constraint W check  (WAVE='R' or WAVE='L'), 
            constraint T check  (TYPE='C' or TYPE='U'),
            constraint M check  (MODE >= 0),
            constraint F check  (FREQ > 0))
            ''')

        self.cursor.execute('''create table PARAMETERS (
            PARAMID      integer primary key autoincrement,
            TYPE         varchar(10) not null,
            LAYER        integer not null,
            UNIT         varchar, 
            constraint U unique (TYPE, LAYER))            
            ''')

        self.cursor.execute('''create table CHAINS (
            CHAINID       integer unique primary key,
            ALGORITHM     varchar default "METROPOLIS")''')

        self.cursor.execute('''create table MODELS (
            MODELID       integer primary key autoincrement,
            CHAINID       int not null references CHAINS(CHAINID),
            WEIGHT        int not null default 1, --weight of the model in the posterior pdf
            NLAYER        int not null, --number of layers including halfspace
            LLK           real not null, --log likelihood of the model
            NITER         int not null default 1) --iteration number for neldermead
            ''')

        self.cursor.execute('''create table PARAMVALUES (
            MODELID       integer references MODELS(MODELID),
            PARAMID       integer references PARAMETERS(PARAMID),
            PARAMVAL      real not null)
            ''')

        self.cursor.execute('''create table DISPVALUES (
            MODELID       integer references MODELS(MODELID),
            POINTID       integer references DISPPOINTS(POINTID),
            DISPVAL       real not null) --null leads to issues when grouping, insert -1.0 instead
            ''')

    def reset(self, nlayer, waves, types, modes, freqs):
        """populate the parameter and disppoints tables"""
        assert np.all([len(waves) == len(_) for _ in (types, modes, freqs)])

        if self.select('select * from PARAMETERS limit 1') is not None:
            raise Exception('PARAMETERS table is not empty')

        if self.select('select * from DISPPOINTS limit 1') is not None:
            raise Exception('DISPPOINTS table is not empty')

        self.begintransaction()
        try:
            for n in range( nlayer):
                if n:
                    self.cursor.execute('insert into PARAMETERS (TYPE, LAYER, UNIT) values (?, ?, ?)', ("Z", n, "KM"))
                self.cursor.execute('insert into PARAMETERS (TYPE, LAYER, UNIT) values (?, ?, ?)', ("VP", n, "KM/S"))
                self.cursor.execute('insert into PARAMETERS (TYPE, LAYER, UNIT) values (?, ?, ?)', ("VS", n, "KM/S"))
                self.cursor.execute('insert into PARAMETERS (TYPE, LAYER, UNIT) values (?, ?, ?)', ("RH", n, "G/CM3"))

            for w, t, m, f in zip(waves, types, modes,  np.round(freqs, 8)):
                self.cursor.execute('''insert into DISPPOINTS (WAVE, TYPE, MODE, FREQ, UNIT) 
                    values (?, ?, ?, ?, "KM/S")''', (w, t, m, f))
            self.commit()
        except:
            self.rollback(crash=True)

    def insert(self, algo, models, datas, weights, llks, parameterizer, datacoder):
        """insert the result of 1 chain (metropolis)"""
        # assume transaction open
        if len(models) == 0:
            warnings.warn('no models provided')
            return

        # give a number to this insertion
        s = self.select('select max(CHAINID) from CHAINS').next()[0]
        if s is None:
            chainid = 0
        else:
            chainid = s + 1

        # determine the parameter ids and point ids
        mdl = models[0]
        ztop, vp, vs, rh = parameterizer.inv(mdl)
        nlayer = len(vs)
        paramids = [] #paramids corresponding to np.concatenate((ztop[1:], vp, vs, rh))
        pointids = [] #pointids corresponding to values

        for i in range(1, nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "Z" and LAYER = ?', (i, )).fetchall()[0][0])
        for i in range(nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "VP" and LAYER = ?', (i, )).fetchall()[0][0])
        for i in range(nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "VS" and LAYER = ?', (i, )).fetchall()[0][0])
        for i in range(nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "RH" and LAYER = ?', (i, )).fetchall()[0][0])

        for w, t, m, f in zip(datacoder.waves, datacoder.types, datacoder.modes,
                              np.round(datacoder.freqs, 8)):
            pointids.append(self.cursor.execute(
                 '''select POINTID from DISPPOINTS where 
                        WAVE=? and TYPE=? and MODE=? and FREQ=?''',
                 (w, t, m, f)).fetchall()[0][0])

        # create the chain
        self.cursor.execute('insert into CHAINS (CHAINID, ALGORITHM) values (?, ?)', (chainid, algo))

        # create the chain
        for niter, (mdl, dat, wgt, llk) in enumerate(zip(models, datas, weights, llks)):
            ztop, vp, vs, rh = parameterizer.inv(mdl)
            paramvals = np.concatenate((ztop[1:], vp, vs, rh))
            dispvals = datacoder.inv(dat)

            dispvals[np.isnan(dispvals)] = -1.0 # convention for nan

            nlayer = len(vs)
            self.cursor.execute('''insert into MODELS (CHAINID, WEIGHT, NLAYER, LLK, NITER)
                values (?, ?, ?, ?, ?)''', (chainid, wgt, nlayer, llk, niter))
            modelid = self.cursor.lastrowid

            self.cursor.executemany('''
                insert into PARAMVALUES (MODELID, PARAMID, PARAMVAL)
                    values(?, ?, ?)
                ''', [(modelid, paramid, paramval)
                      for paramid, paramval in zip(paramids, paramvals)])

            self.cursor.executemany('''
                insert into DISPVALUES (MODELID, POINTID, DISPVAL)
                    values(?, ?, ?)
                ''', [(modelid, pointid, dispval)
                      for pointid, dispval in zip(pointids, dispvals)])

    def get(self, llkmin=None, limit=None, step=None, algo=None):
        """
        get best models from the database, ordered by decreasing likelihood
        input :
            llkmin = None, 0 or float < 0, lowest accepted log-likelihood
            limit = None, 0 or int, max number of distinct models to return
            step = None or int >= 1, yield one model every "step" models
            algo = None or string, select only chains from a given algo ("METROPOLIS" or "NELERMEAD")
                   None means both
        output (generator):
            item = (MODELID, WEIGHT, LLK, NLAYER, (ztop, vp, vs, rh), (waves, types, modes, freqs, values))
                where
                MODELID = int, model number
                WEIGHT = int, weight of the model in the markov chain
                LLK = float, log-likelihood of the model
                NLAYER = int, number of layers in the model
                (ztop, vp, vs, rh) = arrays, depth model
                (waves, types, modes, freqs, values) = arrays with same lengths, dispersion curves
        """

        if llkmin is None: llkmin = 0
        if limit is None: limit = 0
        if step is None: step = 1

        llksql = "where llk > %f" % llkmin if llkmin != 0 else ""
        limitsql = "limit %d" % limit if limit != 0 else ""
        algosql = 'where ALGORITHM = "%s"' % algo if algo is not None else ""

        sql= """
        select MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT,PL, PV, W, T, M, F, DV from
            -- chains, exclude optimization models
                (select CHAINID from CHAINS {algosql}) 

            join        
            -- models
                (select MODELID, CHAINID, WEIGHT, LLK, NLAYER, 
                   group_concat(TYPE) as PT,
                   group_concat(LAYER) as PL, 
                   group_concat(PARAMVAL) as PV
                from 
                    (select * from MODELS {llksql})
                join 
                  (select MODELID, PARAMID, PARAMVAL from PARAMVALUES) using (MODELID)
                join PARAMETERS using (PARAMID)	
                   group by MODELID)
                using (CHAINID)
               
            join
            -- data
                (select MODELID,  
                   group_concat(WAVE) as W, 
                   group_concat(TYPE) as T, 
                   group_concat(MODE) as M, 
                   group_concat(FREQ) as F, 
                   group_concat(DISPVAL) as DV 
                from 
                    (select MODELID from MODELS {llksql})
                join 
                    (select MODELID, POINTID, DISPVAL from DISPVALUES) using (MODELID)
                join DISPPOINTS using (POINTID)	
                group by MODELID)
                using (MODELID)
                        
            order by LLK DESC
            {limitsql}

        """.format(algosql=algosql,
                   llksql=llksql,
                   limitsql=limitsql)

        s = self.select(sql)
        if s is None: return

        for n, (MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT,PL, PV, W, T, M, F, DV) in enumerate(s):
            if n % step: continue

            PT = np.asarray(PT.split(','), '|S2')  # parameter type
            PL = np.asarray(PL.split(','), int)    # parameter layer
            PV = np.asarray(PV.split(','), float)  # parameter value
            Z, VP, VS, RH = [np.zeros(NLAYER, float) for _ in range(4)]  # do not rename variables!!
            for pt, pl, pv in zip(PT, PL, PV):
                #_ = eval(pt)  # Z, VP, VS or RH
                #_[pl] = pv
                eval(pt)[pl] = pv

            W = np.asarray(W.split(','), '|S1')  # waves
            T = np.asarray(T.split(','), '|S1')  # types
            M = np.asarray(M.split(','), int)    # modes
            F = np.asarray(F.split(','), float)  # frequency
            DV = np.asarray(DV.split(','), float)  # dispersion value
            DV[DV == -1.0] = np.nan #convention

            yield MODELID, CHAINID, WEIGHT, LLK, NLAYER, (Z, VP, VS, RH), (W, T, M, F, DV)

    def getzip(self, *args, **kwargs):
        """
        zipped output from self.get
        """
        start = time.time()
        out = list(self.get(*args, **kwargs))
        if not len(out):
            return [np.array([]) for _ in range(5)]

        print ("retrieved %d models in %.6fs " % (len(out), time.time() - start))
        _, chainids, weights, llks, _, ms, ds = zip(*out)

        return np.asarray(chainids, int), \
            np.asarray(weights, float), \
            np.asarray(llks, float), ms, ds

    def getpack(self, *args, **kwargs):
        """
        same as get except that
            - depth model are packed into a depthmodel object
            - dispersion curves are packed into a surf96reader object
        output :
            item = (MODELID, WEIGHT, LLK, NLAYER, dm, sr)
            dm = depthmodel
            sr = surf96reader
        """
        for MODELID, CHAINID, WEIGHT, LLK, NLAYER, (Z, VP, VS, RH), (W, T, M, F, DV) in \
            self.get(*args, **kwargs):

            dm = depthmodel_from_arrays(Z, VP, VS, RH)
            sr = surf96reader_from_arrays(W, T, M, F, DV, None)
            yield MODELID, CHAINID, WEIGHT, LLK, NLAYER, dm, sr

    def getlasts(self, llkmin=None, limit=None, step=None, algo=None):
        """same as self.get to collect the last models generated by the chains (not the best ones)
        """

        if llkmin is None: llkmin = 0
        if limit is None: limit = 0
        if step is None: step = 1

        llksql = "where llk > %f" % llkmin if llkmin != 0 else ""
        limitsql = "limit %d" % limit if limit != 0 else ""
        algosql = 'where ALGORITHM = "%s"' % algo if algo is not None else ""

        sql = """
        select MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT,PL, PV, W, T, M, F, DV from
            -- chains, exclude optimization models
                (select CHAINID from CHAINS {algosql}) 

            join        
            -- models
                (select MODELID, CHAINID, WEIGHT, LLK, NLAYER, 
                   group_concat(TYPE) as PT,
                   group_concat(LAYER) as PL, 
                   group_concat(PARAMVAL) as PV
                from 
                    (select * from MODELS {llksql}
                     order by NITER DESC, CHAINID ASC
                     {limitsql})
                join 
                  (select MODELID, PARAMID, PARAMVAL from PARAMVALUES) using (MODELID)
                join PARAMETERS using (PARAMID)	
                   group by MODELID)
                using (CHAINID)

            join
            -- data
                (select MODELID,  
                   group_concat(WAVE) as W, 
                   group_concat(TYPE) as T, 
                   group_concat(MODE) as M, 
                   group_concat(FREQ) as F, 
                   group_concat(DISPVAL) as DV 
                from 
                    (select MODELID from MODELS {llksql}
                     order by NITER DESC, CHAINID ASC
                     {limitsql})
                join 
                    (select MODELID, POINTID, DISPVAL from DISPVALUES) using (MODELID)
                join DISPPOINTS using (POINTID)	
                group by MODELID)
                using (MODELID)
        order by LLK DESC
        """.format(algosql=algosql,
                   llksql=llksql,
                   limitsql=limitsql)

        s = self.select(sql)
        if s is None: return

        for n, (MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT, PL, PV, W, T, M, F, DV) in enumerate(s):
            if n % step: continue

            PT = np.asarray(PT.split(','), '|S2')  # parameter type
            PL = np.asarray(PL.split(','), int)  # parameter layer
            PV = np.asarray(PV.split(','), float)  # parameter value
            Z, VP, VS, RH = [np.zeros(NLAYER, float) for _ in range(4)]  # do not rename variables!!
            for pt, pl, pv in zip(PT, PL, PV):
                # _ = eval(pt)  # Z, VP, VS or RH
                # _[pl] = pv
                eval(pt)[pl] = pv

            W = np.asarray(W.split(','), '|S1')  # waves
            T = np.asarray(T.split(','), '|S1')  # types
            M = np.asarray(M.split(','), int)  # modes
            F = np.asarray(F.split(','), float)  # frequency
            DV = np.asarray(DV.split(','), float)  # dispersion value
            DV[DV == -1.0] = np.nan  # convention

            yield MODELID, CHAINID, WEIGHT, LLK, NLAYER, (Z, VP, VS, RH), (W, T, M, F, DV)

    def getlastszip(self, *args, **kwargs):
        """
        zipped output from self.getlasts
        """
        start = time.time()
        out = list(self.getlasts(*args, **kwargs))
        if not len(out):
            return [np.array([]) for _ in range(5)]

        print ("retrieved %d models in %.6fs " % (len(out), time.time() - start))
        _, chainids, weights, llks, _, ms, ds = zip(*out)

        return np.asarray(chainids, int), \
               np.asarray(weights, float), \
               np.asarray(llks, float), ms, ds

    def summary(self, head=""):
        """summarize the content of this runfile"""
        s = self.select('''
            select count(*) as NCHAIN, NMODEL, LLKMIN, LLKMAX from 
                CHAINS
                join (select count(*) as NMODEL, MIN(LLK) as LLKMIN, MAX(LLK) as LLKMAX from MODELS)
                 ''')
        if s is None:
            return
        filesize = os.stat(self.sqlitefile).st_size
        NCHAIN, NMODEL, LLKMIN, LLKMAX = s.next()
        if NCHAIN:
            print ("%s%6d chains, %6d models, worst %10f, best %10f, filesize %dM" %
                   (head, NCHAIN, NMODEL, LLKMIN, LLKMAX, filesize / 1024 / 1024))
        else:
            print (head, "no chains")

    def stats(self, head=""):
        """print details about the content of this runfile"""
        s = self.select('''
            select CHAINID, count(*), MAX(LLK) as MLLK from MODELS
                group by CHAINID
                order by MLLK 
            ''')
        if s is None: return
        for CHAINID, N, MLLK in s:
            print ("%schain : %6d  models : %6d maxllk : %f" % (head, CHAINID, N, MLLK))

    def del_sql(self, sql=None):
        """delete models based on a sql command, private"""
        # switch off foerign keys for a moment, otherwise the deletion is fucking slow
        self.cursor.execute('PRAGMA FOREIGN_KEYS = OFF; ')
        try:
            self.cursor.execute('''
                delete from DISPVALUES
                    where MODELID in
                        (select MODELID from MODELS
                            {sql})
                    '''.format(sql=sql))

            self.cursor.execute('''
                delete from PARAMVALUES
                    where MODELID in
                        (select MODELID from MODELS
                            {sql})
                            '''.format(sql=sql))

            self.cursor.execute('''
                delete from MODELS 
                    {sql}
                    '''.format(sql=sql))

        finally:
            self.cursor.execute('''PRAGMA FOREIGN_KEYS = ON; ''')

    def del_bad(self, llkmin):
        sql = "where LLK < {llkmin}".format(llkmin=llkmin)
        self.del_sql(sql=sql)

    def del_chain(self, chainid):
        """delete one or more chains using their chainids,
           deletes only the models, dispvalues and paramvalues,
           preserve entries in the CHAINS table
        """
        if not hasattr(chainid, "__iter__"):
            sql = "where CHAINID = {chainid}".format(chainid=chainid)
        elif len(chainid) == 1:
            sql = "where CHAINID = {chainid}".format(chainid=chainid[0])
        else:
            sql = "where CHAINID in {chainid}".format(chainid=str(tuple(chainid)))
        self.del_sql(sql=sql)


