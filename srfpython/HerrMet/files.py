from tetedenoeud.multipro.multipro8 import Job, MapAsync
from tetedenoeud.utils.asciifile import AsciiFile
from priorpdf import DefaultLogRhoM, LogRhoM_DVS, LogRhoM_DVPDVSDRH, LogRhoM_DVPDVSDRHDPR
from parameterizers import Parameterizer_mZVSPRRH, Parameterizer_mZVSVPRH, Parameterizer_mZVSPRzRHvp, Parameterizer_mZVSPRzRHz
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, depthspace
import numpy as np


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
    print "parameter type : ", p.__class__.__name__
    print "prior type     : ", logRHOM.__class__.__name__
    return p, logRHOM


# -------------------------------------
# run file : crappy, to be optimized
# -------------------------------------
def read_runfile_serial(f):
    with open(f, 'r') as fid:
        nfield = None
        while True:
            l = fid.readline().strip()
            if l == "": break
            if l == "\n" or l.startswith("#"): continue
            l = l.strip('\n').split()
            if nfield is None: nfield = len(l)
            elif not len(l) == nfield:
                print l
                break #line is not full, a run is probably in progress
            chainid, weight, nlayer = np.asarray(l[:3], int)
            llk = float(l[3])
            lmod = np.concatenate(([0.], np.asarray(l[4:4 + 4 * nlayer - 1], float)))

            ldat = l[4 + 4 * nlayer - 1:]
            ndat = len(ldat) / 5

            ztop = lmod[:nlayer]
            vp = lmod[nlayer: 2 * nlayer]
            vs = lmod[2 * nlayer: 3 * nlayer]
            rh = lmod[3 * nlayer: 4 * nlayer]
            waves = ldat[:ndat]
            types = ldat[ndat:2 * ndat]
            modes = np.array(ldat[2 * ndat:3 * ndat], int)
            freqs = np.array(ldat[3 * ndat:4 * ndat], float)
            values = np.array(ldat[4 * ndat:5 * ndat], float)
            yield chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values)


# -------------------------------------
def read_runfile(f, **mapkwargs):
    def gen():
        with open(f, 'r') as fid:
            nfield = None
            while True:
                l = fid.readline().strip()
                if l == "": break
                if l == "\n" or l.startswith("#"): continue
                l = l.strip('\n').split()
                if nfield is None:
                    nfield = len(l)
                elif not len(l) == nfield:
                    print l
                    break  # line is not full, a run is probably in progress
                yield Job(l)

    def fun(l):
        chainid, weight, nlayer = np.asarray(l[:3], int)
        llk = float(l[3])
        lmod = np.concatenate(([0.], np.asarray(l[4:4 + 4 * nlayer - 1], float)))

        ldat = l[4 + 4 * nlayer - 1:]
        ndat = len(ldat) / 5

        ztop = lmod[:nlayer]
        vp = lmod[nlayer: 2 * nlayer]
        vs = lmod[2 * nlayer: 3 * nlayer]
        rh = lmod[3 * nlayer: 4 * nlayer]
        waves = ldat[:ndat]
        types = ldat[ndat:2 * ndat]
        modes = np.array(ldat[2 * ndat:3 * ndat], int)
        freqs = np.array(ldat[3 * ndat:4 * ndat], float)
        values = np.array(ldat[4 * ndat:5 * ndat], float)
        return  chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values)
    with MapAsync(fun, gen(), **mapkwargs) as ma:
        for _, ans, _, _ in ma:
            #print (j[1] - j[0]) / (g[1] - g[0])
            yield ans


# -------------------------------------
def read_runfile_1(f, top=None, topstep=1, **mapkwargs):
    chainids, weights, llks, ms, ds = zip(*list(read_runfile(f, **mapkwargs)))
    chainids, weights, llks = [np.asarray(_) for _ in chainids, weights, llks]
    if top is not None:
        I = np.argsort(llks)[::-1][:top][::topstep]
    else: #means all data
        I = np.argsort(llks)[::-1]
    return chainids[I], weights[I], llks[I], [ms[i] for i in I], [ds[i] for i in I]


# -------------------------------------
def read_runfile_2(f, **mapkwargs):
    with open(f, 'r') as fid:
        nmodels = 1
        while True:
            l = fid.readline()
            if l == "": break
            if l.startswith('#'): continue
            nmodels += 1
    # -------------------
    g = read_runfile(f, **mapkwargs)
    chainid, weight, llk, (ztop, vp, vs, rh), (waves, types, modes, freqs, values) = g.next()
    ZTOP = np.zeros((nmodels, len(ztop)), float) * np.nan
    VP = np.zeros((nmodels, len(vp)), float) * np.nan
    VS   = np.zeros((nmodels, len(vs)), float) * np.nan
    RH = np.zeros((nmodels, len(rh)), float) * np.nan
    WEIGHTS = np.zeros(nmodels, int)
    LLKS = np.zeros(nmodels, int) * np.nan
    # -------------------
    ZTOP[0, :], VP[0, :], VS[0, :], RH[0,:], WEIGHTS[0], LLKS[0] = ztop, vp, vs, rh, weight, llk
    # -------------------
    for n, (chainid, weight, llk, (ztop, vp, vs, rh), _) in enumerate(read_runfile(f, **mapkwargs)):
        ZTOP[n+1, :], \
        VP[n+1, :], \
        VS[n+1, :], \
        RH[n+1, :], \
        WEIGHTS[n+1],\
        LLKS = ztop, vp, vs, rh, weight, llk
    # -------------------
    return ZTOP, VP, VS, RH, WEIGHTS, LLKS


# --------------------
if __name__ == "__main__":
    write_default_paramfile(nlayer=7, zbot=3,
                            type="mZVSPRRH", basedon=None,
                            dvp=None, dvs=None, drh=None, dpr=None)
    A = load_paramfile("_HerrMet.param")
    print A
