from tetedenoeud.database.database import Database
from tetedenoeud.utils.asciifile import AsciiFile
from priorpdf import DefaultLogRhoM, LogRhoM_DVS, LogRhoM_DVPDVSDRH, LogRhoM_DVPDVSDRHDPR
from parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, Parameterizer_mZVSPRzRHvp, Parameterizer_mZVSPRzRHz
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel_from_arrays, depthspace
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
import numpy as np
import time


# -------------------------------------
# param file
# -------------------------------------
def write_default_paramfile(nlayer, zbot, type = "mZVSPRRH", basedon=None, dvp=None, dvs=None, drh=None, dpr=None):
    """create a default parameter file to be customized by user"""
    # ------

    if np.all([_ is None for _ in dvs, dvp, drh, dpr]):
        which = None
    elif dvs is not None and dvp is None and drh is None and dpr is None:
        which = LogRhoM_DVS
    elif dvs is not None and dvp is not None and drh is not None and dpr is None:
        which = LogRhoM_DVPDVSDRH
    elif dvs is not None and dvp is not None and drh is not None and dpr is not None:
        which = LogRhoM_DVPDVSDRHDPR
    else:
        raise NotImplementedError('please specify either dvs alone, or dvp, dvs and drh, or dvp, dvs, drh and dpr')

    # ------
    def write_priortype_header(fid, dvp, dvs, drh):
        if which is None: pass
        elif which is LogRhoM_DVS:
            fid.write('#met PRIORTYPE = "DVS"\n')
            fid.write('#met DVSMIN = %f\n' % dvs[0])
            fid.write('#met DVSMAX = %f\n' % dvs[1])
        elif which is LogRhoM_DVPDVSDRH:
            fid.write('#met PRIORTYPE = "DVPDVSDRH"\n')
            fid.write('#met DVPMIN = %f\n' % dvp[0])
            fid.write('#met DVPMAX = %f\n' % dvp[1])
            fid.write('#met DVSMIN = %f\n' % dvs[0])
            fid.write('#met DVSMAX = %f\n' % dvs[1])
            fid.write('#met DRHMIN = %f\n' % drh[0])
            fid.write('#met DRHMAX = %f\n' % drh[1])
        elif which is LogRhoM_DVPDVSDRHDPR:
            fid.write('#met PRIORTYPE = "DVPDVSDRHDPR"\n')
            fid.write('#met DVPMIN = %f\n' % dvp[0])
            fid.write('#met DVPMAX = %f\n' % dvp[1])
            fid.write('#met DVSMIN = %f\n' % dvs[0])
            fid.write('#met DVSMAX = %f\n' % dvs[1])
            fid.write('#met DRHMIN = %f\n' % drh[0])
            fid.write('#met DRHMAX = %f\n' % drh[1])
            fid.write('#met DPRMIN = %f\n' % dpr[0])
            fid.write('#met DPRMAX = %f\n' % dpr[1])
        else:  raise Exception('programming error')

    # ------
    if type == "mZVSPRRH":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
            prinf = 1.6 * np.ones(nlayer) #r43 * np.ones(nlayer)
            prsup = 2.5 * np.ones(nlayer) #3.5 * np.ones(nlayer)
            rhinf = 2.2 * np.ones(nlayer)
            rhsup = 2.7 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()
            prinf = b.pr().values.copy()
            prsup = b.pr().values.copy()
            rhinf = b.rh.values.copy()
            rhsup = b.rh.values.copy()

        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)] + \
               ["PR%d" % i for i in xrange(nlayer)] + \
               ["RH%d" % i for i in xrange(nlayer)]
        vinfs = np.concatenate((ztopinf, vsinf, prinf, rhinf))
        vsups = np.concatenate((ztopsup, vssup, prsup, rhsup))
        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRRH"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file
    # ----------------------------
    elif type == "mZVSVPRH":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
            vpinf = 0.5 * np.ones(nlayer)
            vpsup = 6.5 * np.ones(nlayer)
            rhinf = 2.2 * np.ones(nlayer)
            rhsup = 2.7 * np.ones(nlayer)

        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()
            vpinf = b.vp.values.copy()
            vpsup = b.vp.values.copy()
            rhinf = b.rh.values.copy()
            rhsup = b.rh.values.copy()


        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)] + \
               ["VP%d" % i for i in xrange(nlayer)] + \
               ["RH%d" % i for i in xrange(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf, vpinf, rhinf))
        vsups = np.concatenate((ztopsup, vssup, vpsup, rhsup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSVPRH"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file
    # ----------------------------
    elif type == "mZVSPRzRHvp":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()

        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRzRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('#met PRz  = "1.0335 * np.exp(-Z / 0.5408) + 1.7310  #some function of Z, Z is in km and growing downward"\n')
            fid.write('#met RHvp = "1.74 * VP ** 0.25 #some function of VP, VP is in km/s, RH is in g/cm3"\n')
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file
    # ----------------------------
    elif type == "mZVSPRzRHz":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            #ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()

        keys = ["-Z%d" % i for i in xrange(1, nlayer)] + \
               ["VS%d" % i for i in xrange(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRzRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('#met PRz  = "1.0335 * np.exp(-Z / 0.5408) + 1.7310  #some function of Z, Z is in km and growing downward"\n')
            fid.write('#met RHz  = "Z * 0. + 2.67 #some function of Z, Z is in km and growing downward"\n')
            write_priortype_header(fid, dvp, dvs, drh)
            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file
    else:
        raise NotImplementedError('no such parameter file type implemented %s' % type)


# -------------------------------------
def load_paramfile(paramfile):
    """initiate one of the parameterizer and prior pdf according to the param file"""
    A = AsciiFile(paramfile)

    # ------------------------
    if A.metadata['TYPE'] == "mZVSVPRH":      p = Parameterizer_mZVSVPRH(A)
    elif A.metadata['TYPE'] == "mZVSPRRH":    p = Parameterizer_mZVSPRRH(A)
    elif A.metadata['TYPE'] == "mZVSPRzRHvp": p = Parameterizer_mZVSPRzRHvp(A)
    elif A.metadata['TYPE'] == "mZVSPRzRHz":  p = Parameterizer_mZVSPRzRHz(A)
    else: raise Exception('could not load %s (TYPE not recognized)' % paramfile)

    # ------------------------
    if not "PRIORTYPE" in A.metadata.keys():
        logRHOM = DefaultLogRhoM(p)
    elif A.metadata['PRIORTYPE'] == "DVS":
        logRHOM = LogRhoM_DVS(p,
                    dvsmin=A.metadata['DVSMIN'],
                    dvsmax=A.metadata['DVSMAX'])
    elif A.metadata['PRIORTYPE'] == "DVPDVSDRH":
        logRHOM = LogRhoM_DVPDVSDRH(p,
                    dvpmin=A.metadata['DVPMIN'],
                    dvpmax=A.metadata['DVPMAX'],
                    dvsmin=A.metadata['DVSMIN'],
                    dvsmax=A.metadata['DVSMAX'],
                    drhmin=A.metadata['DRHMIN'],
                    drhmax=A.metadata['DRHMAX'])
    elif A.metadata['PRIORTYPE'] == "DVPDVSDRHDPR":
        logRHOM = LogRhoM_DVPDVSDRHDPR(p,
                    dvpmin=A.metadata['DVPMIN'],
                    dvpmax=A.metadata['DVPMAX'],
                    dvsmin=A.metadata['DVSMIN'],
                    dvsmax=A.metadata['DVSMAX'],
                    drhmin=A.metadata['DRHMIN'],
                    drhmax=A.metadata['DRHMAX'],
                    dprmin=A.metadata['DPRMIN'],
                    dprmax=A.metadata['DPRMAX'])
    else:
        raise Exception('could not load %s (PRIORTYPE not recognized)' % paramfile)

    # ------------------------
    # print "parameter type : ", p.__class__.__name__
    # print "prior type     : ", logRHOM.__class__.__name__
    return p, logRHOM

# -------------------------------------
# run file : based on sqlite
# -------------------------------------
class RunFile(Database):

    # ----------------------------
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
            CHAINID       integer unique primary key)''')

        self.cursor.execute('''create table MODELS (
            MODELID       integer primary key autoincrement,
            CHAINID       int not null references CHAINS(CHAINID),
            WEIGHT        int not null default 1, --weight of the model in the posterior pdf
            NLAYER        int not null, --number of layers including halfspace
            LLK           real not null ) --log likelyhood of the model
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

    # -------------------
    def reset(self, nlayer, waves, types, modes, freqs):
        """populate the parameter and disppoints tables"""
        assert np.all([len(waves) == len(_) for _ in types, modes, freqs])

        if self.select('select * from PARAMETERS') is not None:
            raise Exception('PARAMETERS table is not empty')

        if self.select('select * from DISPPOINTS') is not None:
            raise Exception('DISPPOINTS table is not empty')

        self.begintransaction()
        try:
            for n in xrange( nlayer):
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

    # -------------------
    def insert(self, models, datas, weights, llks, parameterizer, datacoder):
        """insert the result of 1 chain"""
        # assume transaction open

        # -------------------
        # give a number to this insertion
        s = self.select('select max(CHAINID) from CHAINS').next()[0]
        if s is None:
            chainid = 0
        else:
            chainid = s + 1

        # ------------------
        # determine the parameter ids and point ids
        mdl = models[0]
        ztop, vp, vs, rh = parameterizer.inv(mdl)
        nlayer = len(vs)
        paramids = [] #paramids corresponding to np.concatenate((ztop[1:], vp, vs, rh))
        pointids = [] #pointids corresponding to values

        for i in xrange(1, nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "Z" and LAYER = ?', (i, )).fetchall()[0][0])
        for i in xrange(nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "VP" and LAYER = ?', (i, )).fetchall()[0][0])
        for i in xrange(nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "VS" and LAYER = ?', (i, )).fetchall()[0][0])
        for i in xrange(nlayer):
            paramids.append(self.cursor.execute(
                'select PARAMID from PARAMETERS where TYPE = "RH" and LAYER = ?', (i, )).fetchall()[0][0])

        for w, t, m, f in zip(datacoder.waves, datacoder.types, datacoder.modes,
                              np.round(datacoder.freqs, 8)):
            pointids.append(self.cursor.execute(
                 '''select POINTID from DISPPOINTS where 
                        WAVE=? and TYPE=? and MODE=? and FREQ=?''',
                 (w, t, m, f)).fetchall()[0][0])
        # ------------------
        # create the chain
        self.cursor.execute('insert into CHAINS (CHAINID) values (?)', (chainid, ))

        # create the chain
        for mdl, dat, wgt, llk in zip(models, datas, weights, llks):
            ztop, vp, vs, rh = parameterizer.inv(mdl)
            paramvals = np.concatenate((ztop[1:], vp, vs, rh))
            dispvals = datacoder.inv(dat)


            dispvals[np.isnan(dispvals)] = -1.0 # convention for nan

            nlayer = len(vs)
            self.cursor.execute('''insert into MODELS (CHAINID, WEIGHT, NLAYER, LLK)
                values (?, ?, ?, ?)''', (chainid, wgt, nlayer, llk))
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

    # ----------------------------
    def get(self, llkmin=None, limit=None, step=None):
        """
        get models from the database, ordered by decreasing likelyhood
        input :
            llkmin = None, 0 or float < 0, lowest accepted log-likelihood
            limit = None, 0 or int, max number of distinct models to return
            step = None or int >= 1, yield one model every "step" models
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

        sql= """
        select MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT,PL, PV, W, T, M, F, DV from
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

        """.format(llksql=llksql,
                   limitsql=limitsql)

        # -----------------
        s = self.select(sql)
        if s is None: return

        # -----------------
        for n, (MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT,PL, PV, W, T, M, F, DV) in enumerate(s):
            if n % step: continue

            PT = np.asarray(PT.split(','), '|S2')  # parameter type
            PL = np.asarray(PL.split(','), int)    # parameter layer
            PV = np.asarray(PV.split(','), float)  # parameter value
            Z, VP, VS, RH = [np.zeros(NLAYER, float) for _ in xrange(4)]  # do not rename variables!!
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

    # ----------------------------
    def like_read_run_1(self, llkmin=None, limit=None, step=None):
        """
        zipped output from self.get
        """
        start = time.time()
        out = list(self.get(llkmin=llkmin, limit=limit, step=step))
        print "retrieved %d models in %.6fs " % (len(out), time.time() - start)
        _, chainids, weights, llks, _, ms, ds = zip(*out)

        return np.asarray(chainids, int), \
            np.asarray(weights, float), \
            np.asarray(llks, float), ms, ds

    # ----------------------------
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


# --------------------
if __name__ == "__main__":
    write_default_paramfile(nlayer=7, zbot=3,
                            type="mZVSPRRH", basedon=None,
                            dvp=None, dvs=None, drh=None, dpr=None)
    A = load_paramfile("_HerrMet.param")
    print A



# # -------------------------------------
# # run file : crappy, to be optimized
# # -------------------------------------
# def read_runfile_serial(f):
#     with open(f, 'r') as fid:
#         nfield = None
#         while True:
#             l = fid.readline().strip()
#             if l == "": break
#             if l == "\n" or l.startswith("#"): continue
#             l = l.strip('\n').split()
#             if nfield is None: nfield = len(l)
#             elif not len(l) == nfield:
#                 print l
#                 break #line is not full, a run is probably in progress
#             chainid, weight, nlayer = np.asarray(l[:3], int)
#             llk = float(l[3])
#             lmod = np.concatenate(([0.], np.asarray(l[4:4 + 4 * nlayer - 1], float)))
#
#             ldat = l[4 + 4 * nlayer - 1:]
#             ndat = len(ldat) / 5
#
#             ztop = lmod[:nlayer]
#             vp = lmod[nlayer: 2 * nlayer]
#             vs = lmod[2 * nlayer: 3 * nlayer]
#             rh = lmod[3 * nlayer: 4 * nlayer]
#             waves = ldat[:ndat]
#             types = ldat[ndat:2 * ndat]
#             modes = np.array(ldat[2 * ndat:3 * ndat], int)
#             freqs = np.array(ldat[3 * ndat:4 * ndat], float)
#             values = np.array(ldat[4 * ndat:5 * ndat], float)
#             yield chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values)
#
#
# # -------------------------------------
# def read_runfile(f, **mapkwargs):
#     def gen():
#         with open(f, 'r') as fid:
#             nfield = None
#             while True:
#                 l = fid.readline().strip()
#                 if l == "": break
#                 if l == "\n" or l.startswith("#"): continue
#                 l = l.strip('\n').split()
#                 if nfield is None:
#                     nfield = len(l)
#                 elif not len(l) == nfield:
#                     print l
#                     break  # line is not full, a run is probably in progress
#                 yield Job(l)
#
#     def fun(l):
#         chainid, weight, nlayer = np.asarray(l[:3], int)
#         llk = float(l[3])
#         lmod = np.concatenate(([0.], np.asarray(l[4:4 + 4 * nlayer - 1], float)))
#
#         ldat = l[4 + 4 * nlayer - 1:]
#         ndat = len(ldat) / 5
#
#         ztop = lmod[:nlayer]
#         vp = lmod[nlayer: 2 * nlayer]
#         vs = lmod[2 * nlayer: 3 * nlayer]
#         rh = lmod[3 * nlayer: 4 * nlayer]
#         waves = ldat[:ndat]
#         types = ldat[ndat:2 * ndat]
#         modes = np.array(ldat[2 * ndat:3 * ndat], int)
#         freqs = np.array(ldat[3 * ndat:4 * ndat], float)
#         values = np.array(ldat[4 * ndat:5 * ndat], float)
#         return  chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values)
#     with MapAsync(fun, gen(), **mapkwargs) as ma:
#         for _, ans, _, _ in ma:
#             #print (j[1] - j[0]) / (g[1] - g[0])
#             yield ans
#
#
# # -------------------------------------
# def read_runfile_1(f, top=None, topstep=1, **mapkwargs):
#     chainids, weights, llks, ms, ds = zip(*list(read_runfile(f, **mapkwargs)))
#     chainids, weights, llks = [np.asarray(_) for _ in chainids, weights, llks]
#     if top is not None:
#         I = np.argsort(llks)[::-1][:top][::topstep]
#     else: #means all data
#         I = np.argsort(llks)[::-1]
#     return chainids[I], weights[I], llks[I], [ms[i] for i in I], [ds[i] for i in I]
#
#
# # -------------------------------------
# def read_runfile_2(f, **mapkwargs):
#     with open(f, 'r') as fid:
#         nmodels = 1
#         while True:
#             l = fid.readline()
#             if l == "": break
#             if l.startswith('#'): continue
#             nmodels += 1
#     # -------------------
#     g = read_runfile(f, **mapkwargs)
#     chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values) = g.next()
#     ZTOP = np.zeros((nmodels, len(ztop)), float) * np.nan
#     VP = np.zeros((nmodels, len(vp)), float) * np.nan
#     VS   = np.zeros((nmodels, len(vs)), float) * np.nan
#     RH = np.zeros((nmodels, len(rh)), float) * np.nan
#     WEIGHTS = np.zeros(nmodels, int)
#     LLKS = np.zeros(nmodels, int) * np.nan
#     # -------------------
#     ZTOP[0, :], VP[0, :], VS[0, :], RH[0,:], WEIGHTS[0], LLKS[0] = ztop, vp, vs, rh, weight, llk
#     # -------------------
#     for n, (chainid, weight, llk, (ztop, vp, vs, rh), _) in enumerate(read_runfile(f, **mapkwargs)):
#         ZTOP[n+1, :], \
#         VP[n+1, :], \
#         VS[n+1, :], \
#         RH[n+1, :], \
#         WEIGHTS[n+1],\
#         LLKS = ztop, vp, vs, rh, weight, llk
#     # -------------------
#     return ZTOP, VP, VS, RH, WEIGHTS, LLKS
#
