from srfpython.standalone.database import Database
from srfpython.standalone.asciifile import AsciiFile
from priorpdf import DefaultLogRhoM, LogRhoM_DVS, LogRhoM_DVPDVSDRH, LogRhoM_DVPDVSDRHDPR
from parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, Parameterizer_mZVSPRzRHvp, Parameterizer_mZVSPRzRHz
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthmodel_from_arrays, depthspace
from srfpython.depthdisp.dispcurves import surf96reader_from_arrays
import numpy as np
import time, os


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
            prinf = 1.70 * np.ones(nlayer) #r43 * np.ones(nlayer)
            prsup = 2.15 * np.ones(nlayer) #3.5 * np.ones(nlayer)
            rhinf = 1.8 * np.ones(nlayer)
            rhsup = 3.0 * np.ones(nlayer)
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
            rhinf = 1.8 * np.ones(nlayer)
            rhsup = 3.0 * np.ones(nlayer)

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
            CHAINID       integer unique primary key,
            ALGORITHM     varchar default "METROPOLIS")''')

        self.cursor.execute('''create table MODELS (
            MODELID       integer primary key autoincrement,
            CHAINID       int not null references CHAINS(CHAINID),
            WEIGHT        int not null default 1, --weight of the model in the posterior pdf
            NLAYER        int not null, --number of layers including halfspace
            LLK           real not null, --log likelyhood of the model
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

    # -------------------
    def reset(self, nlayer, waves, types, modes, freqs):
        """populate the parameter and disppoints tables"""
        assert np.all([len(waves) == len(_) for _ in types, modes, freqs])

        if self.select('select * from PARAMETERS limit 1') is not None:
            raise Exception('PARAMETERS table is not empty')

        if self.select('select * from DISPPOINTS limit 1') is not None:
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
    def insert(self, algo, models, datas, weights, llks, parameterizer, datacoder):
        """insert the result of 1 chain (metropolis)"""
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

    # ----------------------------
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
    def getzip(self, *args, **kwargs):
        """
        zipped output from self.get
        """
        start = time.time()
        out = list(self.get(*args, **kwargs))
        if not len(out):
            return [np.array([]) for _ in xrange(5)]

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

    # ----------------------------
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

        # -----------------
        s = self.select(sql)
        if s is None: return

        # -----------------
        for n, (MODELID, CHAINID, WEIGHT, LLK, NLAYER, PT, PL, PV, W, T, M, F, DV) in enumerate(s):
            if n % step: continue

            PT = np.asarray(PT.split(','), '|S2')  # parameter type
            PL = np.asarray(PL.split(','), int)  # parameter layer
            PV = np.asarray(PV.split(','), float)  # parameter value
            Z, VP, VS, RH = [np.zeros(NLAYER, float) for _ in xrange(4)]  # do not rename variables!!
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

    # ----------------------------
    def getlastszip(self, *args, **kwargs):
        """
        zipped output from self.getlasts
        """
        start = time.time()
        out = list(self.getlasts(*args, **kwargs))
        if not len(out):
            return [np.array([]) for _ in xrange(5)]

        print "retrieved %d models in %.6fs " % (len(out), time.time() - start)
        _, chainids, weights, llks, _, ms, ds = zip(*out)

        return np.asarray(chainids, int), \
               np.asarray(weights, float), \
               np.asarray(llks, float), ms, ds

    # ----------------------------
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
            print "%s%6d chains, %6d models, worst %10f, best %10f, filesize %dM" % \
              (head,NCHAIN, NMODEL, LLKMIN, LLKMAX, filesize / 1024 / 1024)
        else:
            print head, "no chains"

    # ----------------------------
    def stats(self, head=""):
        """print details about the content of this runfile"""
        s = self.select('''
            select CHAINID, count(*), MAX(LLK) as MLLK from MODELS
                group by CHAINID
                order by MLLK 
            ''')
        if s is None: return
        for CHAINID, N, MLLK in s:
            print "%schain : %6d  models : %6d maxllk : %f" % (head, CHAINID, N, MLLK)

    # ----------------------------
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
        except:
            raise
        finally:
            self.cursor.execute('''PRAGMA FOREIGN_KEYS = ON; ''')

    # ----------------------------
    def del_bad(self, llkmin):
        sql = "where LLK < {llkmin}".format(llkmin=llkmin)
        self.del_sql(sql=sql)

    # ----------------------------
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


# --------------------
if __name__ == "__main__":
    write_default_paramfile(nlayer=7, zbot=3,
                            type="mZVSPRRH", basedon=None,
                            dvp=None, dvs=None, drh=None, dpr=None)
    A = load_paramfile("_HerrMet.param")
    print A
