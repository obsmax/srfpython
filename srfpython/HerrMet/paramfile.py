from __future__ import print_function
import numpy as np
from srfpython import depthmodel_from_mod96, depthspace, AsciiFile
from srfpython.HerrMet.parameterizers import Parameterizer_mZVSVPRH, Parameterizer_mZVSPRRH, Parameterizer_mZVSPRzRHvp, \
    Parameterizer_mZVSPRzRHz, Parameterizer_mZVSVPvsRHvp
from srfpython.HerrMet.priorpdf import DefaultLogRhoM, LogRhoM_DVS, LogRhoM_DVPDVSDRH, LogRhoM_DVPDVSDRHDPR


"""
"""

# TODO : move this into methods of the parameterizers


def write_default_paramfile(nlayer, zbot,
                            which_parameterizer="mZVSPRRH", basedon=None,
                            dvp=None, dvs=None, drh=None, dpr=None):
    """create a default parameter file to be customized by user"""
    # ------

    if np.all([_ is None for _ in (dvs, dvp, drh, dpr)]):
        which_prior = DefaultLogRhoM
    elif dvs is not None and dvp is None and drh is None and dpr is None:
        which_prior = LogRhoM_DVS
    elif dvs is not None and dvp is not None and drh is not None and dpr is None:
        which_prior = LogRhoM_DVPDVSDRH
    elif dvs is not None and dvp is not None and drh is not None and dpr is not None:
        which_prior = LogRhoM_DVPDVSDRHDPR
    else:
        raise NotImplementedError('please specify either dvs alone, or dvp, dvs and drh, or dvp, dvs, drh and dpr')

    if which_parameterizer == "mZVSPRRH":
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

        keys = ["-Z%d" % i for i in range(1, nlayer)] + \
               ["VS%d" % i for i in range(nlayer)] + \
               ["PR%d" % i for i in range(nlayer)] + \
               ["RH%d" % i for i in range(nlayer)]
        vinfs = np.concatenate((ztopinf, vsinf, prinf, rhinf))
        vsups = np.concatenate((ztopsup, vssup, prsup, rhsup))
        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRRH"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)

            fid.write(which_prior.header(dvp=dvp, dvs=dvs, drh=drh, dpr=dpr))

            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file

    elif which_parameterizer == "mZVSVPRH":
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
            # ztop = np.linspace(0., b.vp.z.max(), nlayer)
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

        keys = ["-Z%d" % i for i in range(1, nlayer)] + \
               ["VS%d" % i for i in range(nlayer)] + \
               ["VP%d" % i for i in range(nlayer)] + \
               ["RH%d" % i for i in range(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf, vpinf, rhinf))
        vsups = np.concatenate((ztopsup, vssup, vpsup, rhsup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSVPRH"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)

            fid.write(which_prior.header(dvp=dvp, dvs=dvs, drh=drh, dpr=dpr))

            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write()  # to screen
        A.write('_HerrMet.param')  # to file

    elif which_parameterizer == "mZVSPRzRHvp":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1] # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
        else:
            b = depthmodel_from_mod96(basedon)
            # ztop = np.linspace(0., b.vp.z.max(), nlayer)
            ztop = depthspace(zbot, nlayer)
            b = b.interp(ztop, interpmethod="weightedstairs")
            ztopinf = -b.vp.z[1:]
            ztopsup = -b.vp.z[1:]
            vsinf = b.vs.values.copy()
            vssup = b.vs.values.copy()

        keys = ["-Z%d" % i for i in range(1, nlayer)] + \
               ["VS%d" % i for i in range(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRzRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('#met PRz  = "def PRz(Z): return 1.0335 * np.exp(-Z / 0.5408) + 1.7310"\n')
            fid.write('#met RHvp = "def RHvp(VP): return 1.74 * VP ** 0.25"\n')

            fid.write(which_prior.header(dvp=dvp, dvs=dvs, drh=drh, dpr=dpr))

            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file

    elif which_parameterizer == "mZVSPRzRHz":
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

        keys = ["-Z%d" % i for i in range(1, nlayer)] + \
               ["VS%d" % i for i in range(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSPRzRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('#met PRz  = "def PRz(Z): return 1.0335 * np.exp(-Z / 0.5408) + 1.7310"\n')
            fid.write('#met RHz  = "def RHz(Z): return Z * 0. + 2.67"\n')

            fid.write(which_prior.header(dvp=dvp, dvs=dvs, drh=drh, dpr=dpr))

            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write() #to screen
        A.write('_HerrMet.param') #to file

    elif which_parameterizer == "mZVSVPvsRHvp":
        if basedon is None:
            ztop = np.linspace(0, zbot, nlayer)
            ztopinf = -ztop[1:]  # deepest side
            ztopsup = -ztop[:-1]  # shallow side
            ztopsup[0] -= 0.001
            vsinf = 0.1 * np.ones(nlayer)
            vssup = 3.5 * np.ones(nlayer)
        else:
            raise NotImplementedError('option -basedon not implemented for mZVSVPvsRHvp')

        keys = ["-Z%d" % i for i in range(1, nlayer)] + \
               ["VS%d" % i for i in range(nlayer)]

        vinfs = np.concatenate((ztopinf, vsinf))
        vsups = np.concatenate((ztopsup, vssup))

        with open('_HerrMet.param', 'w') as fid:
            fid.write('#met TYPE = "mZVSVPvsRHvp"\n')
            fid.write('#met NLAYER = %d\n' % nlayer)
            fid.write('# BROCHER2005\n')
            fid.write(
                '#met VPvs = "def VPvs(VS):'
                'return  0.9409 '
                '+ 2.0947 * VS '
                '- 0.8206 * VS ** 2.0 '
                '+ 0.2683 * VS ** 3.0 '
                '- 0.0251 * VS ** 4.0"  \n')
            fid.write(
                '#met RHvp = "def RHvp(VP):'
                'return  1.6612 * VP '
                '- 0.4721 * VP ** 2.0 '
                '+ 0.0671 * VP ** 3.0 '
                '- 0.0043 * VP ** 4.0 '
                '+ 0.000106 * VP ** 5.0" \n')

            fid.write(which_prior.header(dvp=dvp, dvs=dvs, drh=drh, dpr=dpr))

            fid.write('#fld KEY VINF VSUP\n')
            fid.write('#unt - - -\n')
            fid.write('#fmt %5s %16f %16f\n')
            for k, vinf, vsup in zip(keys, vinfs, vsups):
                fid.write('%s %s %s\n' % (k, vinf, vsup))
        A = AsciiFile('_HerrMet.param')

        A.write()  # to screen
        A.write('_HerrMet.param')  # to file
    else:
        raise NotImplementedError('no such parameter file type implemented %s' % which_parameterizer)


def load_paramfile(paramfile):
    """initiate one of the parameterizer and prior pdf according to the param file"""
    A = AsciiFile(paramfile)

    if A.metadata['TYPE'] == "mZVSVPRH":
        p = Parameterizer_mZVSVPRH(A)
    elif A.metadata['TYPE'] == "mZVSPRRH":
        p = Parameterizer_mZVSPRRH(A)
    elif A.metadata['TYPE'] == "mZVSPRzRHvp":
        p = Parameterizer_mZVSPRzRHvp(A)
    elif A.metadata['TYPE'] == "mZVSPRzRHz":
        p = Parameterizer_mZVSPRzRHz(A)
    elif A.metadata['TYPE'] == "mZVSVPvsRHvp":
        p = Parameterizer_mZVSVPvsRHvp(A)
    else:
        raise Exception('could not load %s (TYPE not recognized)' % paramfile)

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

    print("parameter type : ", p.__class__.__name__)
    print ("prior type     : ", logRHOM.__class__.__name__)
    return p, logRHOM


if __name__ == "__main__":

    write_default_paramfile(nlayer=7, zbot=3,
                            which_parameterizer="mZVSPRRH", basedon=None,
                            dvp=None, dvs=None, drh=None, dpr=None)

    A = load_paramfile("_HerrMet.param")

    print (A)
