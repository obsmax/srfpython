#!/usr/bin/env python
### #!/usr/bin/python2.7

from labex import *
from labex.dispersion.depthmodels import depthmodel, depthmodel1D, depthmodel_from_mod96_txt
# from .
import sys, glob, os, time
import numpy as np
from labex.tools.generic_nobspy import readargv
from labex.tools.stdout import waitbar
from labex.geotools.myproj import Proj
from displininv import displininv, resultdisplay, surf96reader, surf96reader_from_s96_txt, depthmodel_from_mod96_txt, depthpdf
from cube import *

version = "4.1"
default_damp      = 1.0
default_imax      = 10
default_ic        = 0.01
default_sigH      = 0.1
default_sigVS     = 0.1
default_sigPR     = 0.1
default_sigRH     = 0.1
default_smthVS    = 0.0
default_smthPR    = 0.0
default_smthRH    = 0.0
default_nworkers  = 24
#default_prlaw     = "1.0335 * np.exp(-z / 0.5408) + 1.7310" #inferred from GRT1 sonic data
#default_rhlaw     = "1.74 * vp ** 0.25" #Gardner, 1974


help = '''HerrLin V{version}
run linearized depth inversion massively for 3D tomography
-w                change number of workers for parallel computing, default {default_nworkers}
--help, -h        print this message and exits
--target          list of target data (files at surf96 format) to be inverted, 
                  will reproduce file contents into .hlf files with same root names in current directory
                  existing hlf files will be overwriten
    -unc          use it to overwrite uncertainties and set it to constant (in km/s, same for all points in the s96 files)
    -lunc         use it to overwrite uncertainties and set it to constant in log domain (in km/s)
--init            set up initial files for inversion, requires list of hlf files to modify
    -halfspace    change half space parameters in all hlf files (ztop(km), vp(km/s), vs(km/s), rh(km/s)), see also option -lockbot in --inv
    -fwd          if mentionned, I try to forward model the initial model (i.e. compute >forward.000 and >chi2red.000 from >model.000)
                  #initiate with one of the following options : 
                  #1) initialize using a single depthmodel at mod96 format                  
    -m96              depth model to be used as initial model (1 file for all hlf files)
                  #2) initialize with user defined depth model
    -ztop             list of top layer depths, >0, km, first must be 0
    -vp               list of vp in each layer, km/s, same length as ztop
    -vs               list of vs in each layer, km/s, same length as ztop
    -rh               list of rh in each layer, g/cm^3, same length as ztop
                  #3) initialize using one distinct mod96file per hlf files, no more option required
                      mod96files with same root names as hlf files are expected in the current directory
                  #4) initialize using multiple start models for each target
    -mstart           list of multiple start files with same root names as hlffiles 
                      e.g. --target was called to generate 1.hlf and 2.hlf
                      you have too files named 1.foo and 2.bar each containing 100 models to be used as starting models for
                      inversion targets 1 and 2 respectively
                      => will create one hlf file per model and target ./1/start????.hlf and ./2/start????.hlf 
                      use module --inv */start*.hlf to invert each hlf file
                      use module --merge to group solutions corresponding to the same target after inversion
--inv             run iterative linearized inversion (parallelized), requires list hlf files to modify
    -damp         damping, default {default_damp}
    -imax         maximum number of iteration, default {default_imax}
    -ic           interruption condition, default {default_ic}
    -sigH         uncertainty on layer thicknesses (km), default {default_sigH}, 0. means lock layer thicknesses
    -sigVS        uncertainty on Swave velocity (km/s), default {default_sigVS}, 0. means lock VS
    -sigPR        uncertainty on VP/VS ratio (dimentionless), default {default_sigPR}, 0. means lock VP/VS
    -sigRH        uncertainty on RH ratio (g/cm^3), default {default_sigRH}, 0. means lock density
    -smthVS       smoothing distance in layer number, default {default_smthVS}
    -smthPR       smoothing distance in layer number, default {default_smthPR}
    -smthRH       smoothing distance in layer number, default {default_smthRH}
    -lockbot      if mentionned, try to lock halfspace vp, vs and rh by increasing damping for these layers 
                  cannot lock the depth of the half space, default False
--smooth          increase the number of layers in last model, recompute data
                  run it after --init or --inv
                  requires list of hlf files to modify
--rollback        remove all iterations (keep model.000, forward.000, chi2red.000 if exists)
                  requires list of hlf files to modify
--count           count the number of iterations found in hlf files, * means that all forwards have been computed
                  requires list of hlf files to modify
--merge           merge solutions for multi start inversion 
                  take the median of models rootname/start*.hlf
                  and write it as rootname/median.mod96 rootname/p16.mod96 rootname/p84.mod96
                  requires list of hlf files to read
--show            display hlf files and pause
                  requires list of hlf files to modify
    -plim         set period display range
    -vlim         set velocity display range
    -zlim         set depth display range
    -vslim        set vs display range
    -vplim        set vp display range
    -prlim        set vp/vs display range
    -rhlim        set rh display range
    -nstep        display only several iteration over nstep 
                  0 means start only, 
                 -1 means last only, 
                 -2 means 0 & -1, 
                  1 means all (default), 
                  >1 means 0 and -1 and 1 model over nstep
--savefig         generate figures and save it at png format (parallelized)
    same as --show, create png files instead of showing figures
--cube            store inverted files into a 3D model
                  requires list of hlf or m96 files to read
                  requires file names to include _lon* and _lat*
    -dpt          depth array, zmin, zmax, zrange (>0, km)
    -utm          utm zone, e.g. 33N
    -isec         interactive sections, specify depths at which to compute maps
#--------------------------------------
example : 
HerrLin4_1 --target /path/to/*.surf96 -unc 0.1
HerrLin4_1 --init ./*.hlf -m96 /path/to/model0.mod -halfspace 3.0 5.8 3.35 2.67 -fwd
HerrLin4_1 -w 24 --inv ./*.hlf -damp 1.0 -sigH 0.1 -sigVS 0.1 -sigPR 0.01 -sigRH 0.01 -smthVS 0.0 -smthPR 0.0 -smthRH 0.0 -imax 10 -ic 0.01 -lockbot
HerrLin4_1 -w 24 --smooth *.hlf
HerrLin4_1 -w 24 --savefig *partof*.hlf -plim 0.5 5. -vlim 0.4 3. -zlim 0. 3.1 -vslim 0.5 4.0 -vplim 1.0 7.0 -prlim 1.2 3.0 -rhlim 1.8 3.0 
HerrLin4_1 -w 24 --cube *.hlf -dpt 0. 3.1 0.05 -utm 32N -isec 0.1 1. 2.'''.format(\
    version           = version,
    default_nworkers  = default_nworkers,
    default_damp      = default_damp,
    default_imax      = default_imax,
    default_ic        = default_ic,
    default_sigH      = default_sigH,
    default_sigVS     = default_sigVS,
    default_sigPR     = default_sigPR,
    default_sigRH     = default_sigRH,
    default_smthVS    = default_smthVS,
    default_smthPR    = default_smthPR,
    default_smthRH    = default_smthRH)

    
#    default_prlaw     = default_prlaw,
#    default_rhlaw     = default_rhlaw)
#---------------------------------------------------
def strinfile(f, s):
    with open(f, 'r') as fid:
        while True:
            l = fid.readline()
            if l == "":  return False
            elif s in l: return True
#---------------------------------------------------
def s962str(s96):
    if not os.path.exists(s96): raise Exception('%s not found' % s96)
    assert s96.endswith('.surf') or s96.endswith('.surf96')
    L = []
    with open(s96, 'r') as fid: 
        for l in fid.readlines():
            l = l.strip("\n").strip()
            if l.startswith('SURF96'): L.append(l)
            else: print "in %s : line %s not understood" % (s96, l)
    return "\n".join(L)
#---------------------------------------------------
def m962str(m96):
    if not os.path.exists(m96): raise Exception('%s not found' % m96)
    assert m96.endswith('.mod') or m96.endswith('.mod96')
    L = []
    with open(m96, 'r') as fid: 
        for l in fid.readlines():
            l = l.strip("\n").rstrip()
            if l != "": L.append(l)
    return "\n".join(L)
#---------------------------------------------------
def readfloatafter(s, key):
    #read float values from a string after keyword key
    if not key in s:
        raise Exception('keyword %s not found in string %s' % (key, s))
    s = s.split(key)[-1]
    while len(s):
        try:     return float(s)
        except:  s = s[:-1]
    raise Exception('error')
#---------------------------------------------------
def readmstart(filename):
    """read multiple start file
    yields all depth model at mod96 format found in a single ascii file
    starting signal = MODEL.01
    ending   signal = layer with thickness 0
    """
    with open(filename, 'r') as fid:
        l = fid.readline()
        while True:
            if l == "" : break
            l = l.strip('\n').strip()
            if l.startswith('MODEL.01'):
                m = [l]
                while True:
                    l = fid.readline()
                    if l == "" : break
                    l = l.strip('\n').strip()
                    m.append(l)                    
                    if len(l.split()) == 10:
                        try: 
                            H = float(l.split()[0])
                            if H == 0. : break
                        except KeyboardInterrupt: raise
                        except: pass
                yield "\n".join(m)
            l = fid.readline()


#---------------------------------------------------
class HerrLinFile(object):
    """
    manage hlf files, an ascii files with inversion history

    >target surf96filename
    data to invert at surf96 format
    >model.000 mod96filename
    initial model at mod96 format
    >forward.000
    dispersion data forward modeled from model.000 at surf96 format
    >chi2red.000
    reduced data misfit corresponding to model.000
    >model.001
    ...
    """
    #-----------------------------------------------
    def __init__(self, hlf):
        assert isinstance(hlf, str)
        assert hlf.endswith('.hlf')
        self.hlf = hlf
        self.target = None #a tuple like (targetfilename, s96str)
        self.init   = None #name of initfile
        self.models   = [] #list of models (m96str), first = init
        self.forwards = [] #list of dispersion curves (s96str)
        self.chi2reds = [] #list of floats
        if os.path.exists(self.hlf):
            with open(self.hlf, 'r') as fid:
                l = fid.readline()
                while True:
                    if l == "": break
                    l = l.strip().rstrip('\n')
                    if l.startswith('>'):
                        if ">target" in l:
                            targetfile = l.split('>target')[-1].strip() 
                            L = []
                            while True:
                                l = fid.readline()
                                if l == "" or l.startswith(">"): break
                                L.append(l)
                            self.target = (targetfile, "".join(L).strip('\n'))
                            continue #next line already read in
                        elif ">model" in l:
                            nmodel = int(l.split('>model.')[-1].split()[0])
                            if not nmodel:  self.init = l.split()[-1]
                            assert len(self.models) == nmodel
                            L = []
                            while True:
                                l = fid.readline()
                                if l == "" or l.startswith(">"): break
                                L.append(l)
                            self.models.append("".join(L).strip('\n'))
                            continue #next line already read in
                        elif ">forward" in l:
                            nmodel = int(l.split('>forward.')[-1].split()[0])
                            assert len(self.forwards) == nmodel
                            L = []
                            while True:
                                l = fid.readline()
                                if l == "" or l.startswith(">"): break
                                L.append(l)
                            self.forwards.append("".join(L).strip('\n'))
                            continue #next line already read in
                        elif ">chi2red" in l:
                            nmodel = int(l.split('>chi2red.')[-1].split()[0])
                            chi2red = float(fid.readline().strip().strip('\n'))
                            self.chi2reds.append(chi2red)
                    l = fid.readline()
            assert len(self.forwards) == len(self.chi2reds)
            assert len(self.forwards) == len(self.models) or \
                   len(self.forwards) + 1 == len(self.models)

#        print self.chi2reds
#        for m, f, c in zip(self.models, self.forwards, self.chi2reds):
#            print m
#            print f
#            print c
    #-----------------------------------------------
    def __str__(self):
        s = ""
        if self.target is None: return s
        else:
            targetfile, L = self.target
            s += ">target %s\n%s\n" % (targetfile, L)
            if not len(self.models): return s
            else:
                for nmodel in xrange(len(self.models)):
                    s += ">model.%03d %s\n%s\n" % (nmodel, "" if nmodel else self.init, self.models[nmodel])
                    if len(self.forwards) > nmodel:
                        s += ">forward.%03d\n%s\n" % (nmodel, self.forwards[nmodel])
                        s += ">chi2red.%03d\n%f\n" % (nmodel, self.chi2reds[nmodel])               
        return s#.strip('\n').strip()
    #-----------------------------------------------
    def set_target(self, filename, data):
        if self.target is not None:
            raise Exception('%s has already a target' % self.hlf)

        self.target = (filename, data)
        with open(self.hlf, 'a') as fout:
            fout.write('>target %s\n' % filename)
            fout.write('%s\n' % data)
    #-----------------------------------------------
    def set_init(self, m96):
        if self.target is None:
            raise Exception('%s has no target, please set it first' % self.hlf)
        if not len(self.models) == 0:
            raise Exception('%s has already an initial model' % self.hlf)

        if os.path.exists(m96): #m96 is a file name
            L = m962str(m96)
            self.init = m96
        else: #m96 is a file content
            L = m96
            self.init = "_"    

        self.models.append(L)

        with open(self.hlf, 'a') as fout:
            fout.write('>model.000 %s\n' % self.init)
            fout.write('%s\n' % L)
    #---------------------------------------------------
    def change_init_halfspace(self, ZTOP, VP, VS, RH):
        """change half space parameters in init modeled 
           must be called after set_init and before running inversion 
        """
        if self.target is None:
            raise Exception('%s has no target, please set it first' % self.hlf)
        if len(self.models) == 0:
            raise Exception('%s has already no initial model, please set it first' % self.hlf)
        elif len(self.models) > 1:
            raise Exception('%s inversion started already, please roll it back first' % self.hlf)
        assert not hasattr(ZTOP, "__iter__")
        assert not hasattr(VP,   "__iter__")
        assert not hasattr(VS,   "__iter__")
        assert not hasattr(RH,   "__iter__")
        
        dm = depthmodel_from_mod96_txt(self.models[0])
        ztop, vp, vs, rh = dm.vp.z, dm.vp.values, dm.vs.values, dm.rh.values
        if ZTOP > ztop[-1]:
            ztop = np.concatenate((ztop, [ZTOP]))
            vp   = np.concatenate((vp,   [VP]))
            vs   = np.concatenate((vs,   [VS]))
            rh   = np.concatenate((rh,   [RH]))
        elif ZTOP == ztop[-1]:
            vp[-1] = VP
            vs[-1] = VS
            rh[-1] = RH
        else:
            i = np.searchsorted(ztop, ZTOP)
            ztop = np.concatenate((ztop[:i], [ZTOP]))
            vp   = np.concatenate((vp[:i],   [VP]))
            vs   = np.concatenate((vs[:i],   [VS]))
            rh   = np.concatenate((rh[:i],   [RH]))
        dm = depthmodel(\
            depthmodel1D(ztop, vp),
            depthmodel1D(ztop, vs),
            depthmodel1D(ztop, rh))
        self.models[0] = str(dm)
        self.init += "*"

        with open(self.hlf, 'w') as fout:
            fout.write(str(self))#.rstrip().rstrip('\n'))

    #-----------------------------------------------
    def addmodel(self, ZTOP, VP, VS, RH):
        assert len(self.chi2reds) == len(self.forwards) == len(self.models)
        dm = depthmodel(\
            depthmodel1D(ZTOP, VP),
            depthmodel1D(ZTOP, VS),
            depthmodel1D(ZTOP, RH))
        self.models.append(dm.__str__().strip('\n').strip())
    #-----------------------------------------------
    def addforward(self, chi2red, waves, types, modes, freqs, values, dvalues = None):
        assert len(self.chi2reds) == len(self.forwards) == len(self.models) - 1
        if dvalues is None:
            L = []
            for w, t, m, f, v in zip(waves, types, modes, freqs, values):
                L.append('SURF96 %s %s X %d %f %f %f' % (w, t, m, 1. / f, v, 0.1))
            L = "\n".join(L)
        else:
            raise NotImplementedError('')
        self.chi2reds.append(chi2red)
        self.forwards.append(L)
    #-----------------------------------------------         
    def show(self, rd = None, plim = None, vlim = None, zlim = None, vplim = None, vslim = None, rhlim = None, prlim = None, nstep = 1, cmap = tej()):
        target = surf96reader_from_s96_txt(self.target[1])

        if rd is None:rd = resultdisplay()
        rd.plotdisp(color =  "k", linewidth = 3, alpha = 1., *target.wtmfvd())       
        
        N = len(self.forwards)
        clrs = [value2color(n /float(N), cmap = cmap) for n in np.arange(N)]

        def subshow(i, linewidth = 1, alpha = 1, linestyle = "-"):
            model = depthmodel_from_mod96_txt(self.models[i])
            model = (model.vp.z, model.vp.values, model.vs.values, model.rh.values)
            data  = surf96reader_from_s96_txt(self.forwards[i])
            data  = data.wtmfv()

            rd.plotmodel(color = clrs[i], linewidth = linewidth, alpha = alpha, linestyle = linestyle, *model)
            rd.plotdisp(color =  clrs[i], linewidth = linewidth, alpha = alpha, linestyle = linestyle, *data)
            rd.axconv.plot(np.arange(N)[i], self.chi2reds[i], 'o', color = clrs[i])

        if nstep == -1: #last model only
            subshow(-1, linewidth = 1, alpha = 0.2)
        elif nstep == 0: #first model only
            subshow(0, linewidth = 1, alpha = 0.2)
        elif nstep == -2: #first and last models only
            subshow(0, linewidth = 1, alpha = 0.2)
            subshow(-1, linewidth = 1, alpha = 0.2)
        else:
            subshow(0, linewidth = 1, alpha = 1)
            #display intermediate steps
            for i in np.arange(0, N, nstep)[1:-1]:
                subshow(i, linewidth = 1, alpha = 0.4)
            subshow(-1, linewidth = 3, alpha = 1)


        if plim is not None: rd.set_plim(plim)
        if vlim is not None: rd.set_vlim(vlim)
        if zlim is not None: rd.set_zlim(zlim)

        if vplim is not None: rd.axvp.set_xlim(vplim)
        if vslim is not None: rd.axvs.set_xlim(vslim)
        if prlim is not None: rd.axpr.set_xlim(prlim)
        if rhlim is not None: rd.axrh.set_xlim(rhlim)
    
        rd.grid()
        rd.tick()
        #rd.fig.suptitle(self.target[0])
        rd.fig.suptitle(self.hlf)

#---------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print help
        sys.exit()
    argv = readargv()
    nworkers = int(argv['w'][0]) if "w" in argv.keys() else default_nworkers

    #-------------------------------------
    if "h" in argv.keys() or "help" in argv.keys():
        print help
        sys.exit()
    #-------------------------------------
    elif "v" in argv.keys() or "version" in argv.keys():
        print "version : %s" % version
        sys.exit()
    #-------------------------------------
    elif "target" in argv.keys():
        #s96s = argv['s96']
        s96s = argv['target']
        hlfs = []
        for s96 in s96s:
            assert s96.endswith(".surf96")
            assert os.path.exists(s96)
            hlf = s96.split('/')[-1].replace('.surf96', '.hlf')
            hlfs.append(hlf)
        assert len(hlfs) == len(np.unique(hlfs)) #hlf files must have distinct names

        for s96, hlf in zip(s96s, hlfs):
            if os.path.exists(hlf):
                os.system('rm -f %s' % hlf)
            print "%s > %s" % (s96, hlf) 

            #-----------------
            s = surf96reader(s96)
            if "unc" in argv.keys():
                s.data['DVALUE'] = float(argv['unc'][0]) * np.ones_like(s.data['VALUE']) #set constant uncertainty
            elif "lunc" in argv.keys():
                s.data['DVALUE'] = float(argv['lunc'][0]) * s.data.VALUE #set constant uncertainty to the logarithm of the value
            
            #-----------------
            h = HerrLinFile(hlf)
            h.set_target(filename = s96, data = s.__str__())
        sys.exit()
    #-------------------------------------
    elif "init" in argv.keys():
        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['init']

        if "m96" in argv.keys():
            """initiate all inversions from given m96 file"""
            assert len(argv['m96'])  == 1
            for k in "ztop vp vs rh".split() :
                assert not k in argv.keys()

            m96 = argv['m96'][0]
            assert os.path.exists(m96)
            assert m96.endswith('.mod96') or m96.endswith('.mod')

            for hlf in hlfs:
                h = HerrLinFile(hlf)
                print "%s >> %s" % (m96, hlf) 
                h.set_init(m96)#str)
                if "halfspace" in argv.keys():
                    ztop, vp, vs, rh = argv['halfspace']
                    h.change_init_halfspace(ztop, vp, vs, rh)

        elif "ztop" in argv.keys():
            """initiate all inversions using ztop, vs, vp, rh"""
            ztop = argv['ztop']
            vs   = argv['vs']
            vp   = argv['vp']
            rh   = argv['rh']
            assert len(ztop) == len(vs) == len(vp) == len(vs)

            dm = depthmodel(vp = depthmodel1D(ztop, vp),
                            vs = depthmodel1D(ztop, vs),
                            rh = depthmodel1D(ztop, rh))
            
            m96str = str(dm)
            for hlf in hlfs:
                h = HerrLinFile(hlf)
                print "ztop, vs, vp, rh >> %s" % hlf
                h.set_init(m96str.strip().strip('\n').strip())
                if "halfspace" in argv.keys():
                    ztop, vp, vs, rh = argv['halfspace']
                    h.change_init_halfspace(ztop, vp, vs, rh)

        elif "mstart" in argv.keys():
            """target has been called with 
            1.surf96 => 1.hlf
            2.surf96 => 2.hlf
            
            call init with -mstart 1.foo 2.bar
            all models found in 1.foo will be attributed to 1.hlf as ./1/start0000.hlf, ./1/start0001.hlf, ...
            all models found in 2.bar will be attributed to 2.hlf as ./2/start*.hlf
            """
            newhlfs = []
            for hlf in hlfs:
                hlfrootname = hlf.split('/')[-1].split('.hlf')[0]
                os.system('mkdir %s' % hlfrootname)
                print hlf
                for mstart in argv['mstart']:
                    mstartrootname = ".".join(mstart.split('/')[-1].split('.')[:-1])
                    if mstartrootname == hlfrootname:
                        print "    ", mstart
                        for n, mstart in enumerate(readmstart(mstart)):
                            newhlf = "%s/start%04d.hlf" % (hlfrootname, n)
                            os.system('cp %s %s' % (hlf, newhlf))
                            newh = HerrLinFile(newhlf)
                            newh.set_init(mstart)
                            if "halfspace" in argv.keys():
                                ztop, vp, vs, rh = argv['halfspace']
                                newh.change_init_halfspace(ztop, vp, vs, rh)
                            newhlfs.append(newhlf)

            hlfs = newhlfs #needed by option -fwd below
            
        else:
            """find a mod96 file for each hlf file and use it as a starting model"""
            for k in "ztop vp vs rh".split() :
                assert not k in argv.keys()

            m96s = []            
            for hlf in hlfs:
                m96 = hlf.replace('.hlf', '.mod96')
                if not os.path.exists(m96):
                    raise Exception('%s not found (use -m96 to use similar starting model for all hlf files)' % m96)
                m96s.append(m96)
            for hlf, m96 in zip(hlfs, m96s):
                h = HerrLinFile(hlf)
                print "%s >> %s" % (m96, hlf) 

                h.set_init(m96)
                if "halfspace" in argv.keys():
                    ztop, vp, vs, rh = argv['halfspace']
                    h.change_init_halfspace(ztop, vp, vs, rh)

        #-------------------------- 
        if "fwd" in argv.keys():
            def fun(h, *args, **kwargs):
                return h, list(displininv(*args, **kwargs))

            def gen():
                for hlf in hlfs:
                    h = HerrLinFile(hlf)
                    assert len(h.models) >= 1

                    _, target = h.target
                    nlast, last = len(h.models) - 1, h.models[-1]

                    yield Job(h = h, \
                        s96               = target, 
                        mod96             = last, 
                        damping           = 1.,
                        imax              = 0,
                        interruptconditon = 0.1,
                        sigH              = 0.1,
                        sigVS             = 0.1,
                        sigPR             = 0.1, 
                        sigRH             = 0.1,
                        smoothVS          = 0.,
                        smoothPR          = 0.,
                        smoothRH          = 0.)

            wb = waitbar('progress (forward)')
            with MapAsync(fun, gen(), Nworkers = nworkers, RaiseIfError = True) as ma:
                for nnn, (jobid, (h, o), _, _) in enumerate(ma):

                    for n, (chi2red, (ZTOP, VP, VS, RH), (waves, types, modes, freqs, values)) in enumerate(o):
                        if len(h.forwards) == len(h.models):
                            h.addmodel(ZTOP, VP, VS, RH)
                            h.addforward(chi2red, waves, types, modes, freqs, values)
                        else:
                            h.addforward(chi2red, waves, types, modes, freqs, values)
                    #print "+%d iterations >> %s" % (n, h.hlf)
                    with open(h.hlf, 'w') as fid:
                        fid.write(str(h))
                    wb.refresh(nnn / float(len(hlfs)))
                wb.close()


        sys.exit()

    #-------------------------------------
    elif "smooth" in argv.keys():
        hlfs = argv['smooth']
        nsmooth = 1
        def gen():
            for hlf in hlfs:
                h = HerrLinFile(hlf)
                if not len(h.models) >= 1:
                    raise Exception('%s %d' % (h.hlf, len(h.models)))
                yield Job(h)
        def fun(h):
            dm = depthmodel_from_mod96_txt(h.models[-1])
            for _ in xrange(nsmooth):
                pr = dm.pr(); pr.smooth()
                dm.vs.smooth()
                dm.rh.smooth()
                dm = depthmodel_from_vs_pr_rh(dm.vs, pr, dm.rh)
            h.models.append(str(dm))
            return h

        wb = waitbar('progress (smooth)', reevaluatespeed = 3600.)
        with MapAsync(fun, gen(), Nworkers = nworkers, RaiseIfError = True) as ma:
            for nnn, (jobid, h, _, _) in enumerate(ma):
                with open(h.hlf, 'w') as fid:
                    fid.write(str(h))
                wb.refresh(nnn / float(len(hlfs)))
            wb.close()

        #-------------------------- 
        if True: #"fwd" in argv.keys():
            def fun(h, *args, **kwargs):
                return h, list(displininv(*args, **kwargs))

            def gen():
                for hlf in hlfs:
                    h = HerrLinFile(hlf)
                    assert len(h.models) >= 1

                    _, target = h.target
                    nlast, last = len(h.models) - 1, h.models[-1]

                    yield Job(h = h, \
                        s96               = target, 
                        mod96             = last, 
                        damping           = 1.,
                        imax              = 0,
                        interruptconditon = 0.1,
                        sigH              = 0.1,
                        sigVS             = 0.1,
                        sigPR             = 0.1, 
                        sigRH             = 0.1,
                        smoothVS          = 0.,
                        smoothPR          = 0.,
                        smoothRH          = 0.)

            wb = waitbar('progress (forward)')
            with MapAsync(fun, gen(), Nworkers = nworkers, RaiseIfError = True) as ma:
                for nnn, (jobid, (h, o), _, _) in enumerate(ma):

                    for n, (chi2red, (ZTOP, VP, VS, RH), (waves, types, modes, freqs, values)) in enumerate(o):
                        if len(h.forwards) == len(h.models):
                            h.addmodel(ZTOP, VP, VS, RH)
                            h.addforward(chi2red, waves, types, modes, freqs, values)
                        else:
                            h.addforward(chi2red, waves, types, modes, freqs, values)
                    #print "+%d iterations >> %s" % (n, h.hlf)
                    with open(h.hlf, 'w') as fid:
                        fid.write(str(h))
                    wb.refresh(nnn / float(len(hlfs)))
                wb.close()

        sys.exit()

    #-------------------------------------
    elif "inv" in argv.keys():
        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['inv']

        print '''invert with :
        damping           = ''' + str(argv['damp'][0]         if "damp"  in argv.keys() else default_damp)   + '''
        imax              = ''' + str(int(argv['imax'][0])    if "imax"  in argv.keys() else default_imax)   + '''
        interruptconditon = ''' + str(argv['ic'][0]           if "ic"    in argv.keys() else default_ic)     + '''
        sigH              = ''' + str(argv['sigH'][0]         if "sigH"  in argv.keys() else default_sigH)   + '''
        sigVS             = ''' + str(argv['sigVS'][0]        if "sigVS" in argv.keys() else default_sigVS)  + '''
        sigPR             = ''' + str(argv['sigPR'][0]        if "sigPR" in argv.keys() else default_sigPR)  + '''
        sigRH             = ''' + str(argv['sigRH'][0]        if "sigRH" in argv.keys() else default_sigRH)  + '''
        smoothVS          = ''' + str(argv['smthVS'][0]       if "smthVS"  in argv.keys() else default_smthVS)   + '''
        smoothPR          = ''' + str(argv['smthPR'][0]       if "smthPR"  in argv.keys() else default_smthPR)   + '''
        smoothRH          = ''' + str(argv['smthRH'][0]       if "smthRH"  in argv.keys() else default_smthRH)      + '''
        lockbot           = ''' + str("lockbot" in argv.keys())

        def fun(h, *args, **kwargs):
            return h, list(displininv(*args, **kwargs))

        def gen():
            for hlf in hlfs:
                h = HerrLinFile(hlf)
                if not len(h.models) >= 1:
                    raise Exception('%s %d' % (h.hlf, len(h.models)))

                _, target = h.target
                nlast, last = len(h.models) - 1, h.models[-1]

                yield Job(h = h, \
                    s96               = target, 
                    mod96             = last, 
                    damping           = argv['damp'][0]       if "damp"    in argv.keys() else default_damp,
                    imax              = int(argv['imax'][0])  if "imax"    in argv.keys() else default_imax,
                    interruptconditon = argv['ic'][0]         if "ic"      in argv.keys() else default_ic,
                    sigH              = argv['sigH'][0]       if "sigH"    in argv.keys() else default_sigH,
                    sigVS             = argv['sigVS'][0]      if "sigVS"   in argv.keys() else default_sigVS,
                    sigPR             = argv['sigPR'][0]      if "sigPR"   in argv.keys() else default_sigPR,
                    sigRH             = argv['sigRH'][0]      if "sigRH"   in argv.keys() else default_sigRH,
                    smoothVS          = argv['smthVS'][0]     if "smthVS"  in argv.keys() else default_smthVS,
                    smoothPR          = argv['smthPR'][0]     if "smthPR"  in argv.keys() else default_smthPR,
                    smoothRH          = argv['smthRH'][0]     if "smthRH"  in argv.keys() else default_smthRH,
                    lockbot           = 'lockbot' in argv.keys())
#                    prlaw             = str(argv['prlaw'][0]) if "prlaw" in argv.keys() else default_prlaw,
#                    rhlaw             = str(argv['rhlaw'][0]) if "rhlaw" in argv.keys() else default_rhlaw)

        wb = waitbar('progress (inversion)', reevaluatespeed = 3600.)
        with MapAsync(fun, gen(), Nworkers = nworkers, RaiseIfError = True) as ma:
            for nnn, (jobid, (h, o), _, _) in enumerate(ma):
                firstdone =  len(h.forwards) == len(h.models)
                for n, (chi2red, (ZTOP, VP, VS, RH), (waves, types, modes, freqs, values)) in enumerate(o):
                    if not n and firstdone: continue
                    if len(h.forwards) == len(h.models):
                        h.addmodel(ZTOP, VP, VS, RH)
                        h.addforward(chi2red, waves, types, modes, freqs, values)
                    else:
                        h.addforward(chi2red, waves, types, modes, freqs, values)
                #print "+%d iterations >> %s" % (n, h.hlf)
                with open(h.hlf, 'w') as fid:
                    fid.write(str(h))
                wb.refresh(nnn / float(len(hlfs)))
            wb.close()
        sys.exit()
    #-------------------------------------
    elif "rollback" in argv.keys():
        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['rollback']
        for hlf in hlfs:
            h = HerrLinFile(hlf)
            if len(h.models):
                h.models = [h.models[0]]
                if len(h.forwards):
                    h.forwards = [h.forwards[0]]
                    h.chi2reds = [h.chi2reds[0]]
                with open(h.hlf, 'w') as fid:
                    fid.write(str(h))                
        sys.exit()
    #-------------------------------------
    elif "count" in argv.keys():
        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['count']
        for hlf in hlfs:
            h = HerrLinFile(hlf)
            if True:#len(h.models) > 1:
                print "%50s %6s %10f" % (hlf, str(len(h.models)) + ("*" if len(h.forwards) == len(h.models) else ""), h.chi2reds[-1] if len(h.chi2reds) else np.nan)
        sys.exit()
    #-------------------------------------
    elif "show" in argv.keys():
        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['show']
        rd = resultdisplay(fig = None)
        if "m96" in argv.keys():
            dms = [depthmodel_from_mod96(f) for f in argv['m96']]
        else: dms = []
                
        for hlf in hlfs:
            if not "holdon" in argv.keys(): rd.cla()
            h = HerrLinFile(hlf)
            h.show(rd = rd, 
                plim  = argv['plim']  if 'plim'  in argv.keys() else None,
                vlim  = argv['vlim']  if 'vlim'  in argv.keys() else None,
                zlim  = argv['zlim']  if 'zlim'  in argv.keys() else None,
                vplim = argv['vplim'] if 'vplim' in argv.keys() else None,
                vslim = argv['vslim'] if 'vslim' in argv.keys() else None,
                rhlim = argv['rhlim'] if 'rhlim' in argv.keys() else None,
                prlim = argv['prlim'] if 'prlim' in argv.keys() else None,
                nstep = int(argv['nstep'][0]) if 'nstep' in argv.keys() else 1)

            for dm in dms: 
                rd.plotmodel(dm.vp.z, dm.vp.values, dm.vs.values, dm.rh.values, color ="k", linewidth = 3, alpha = 1.)

            if "holdon" in argv.keys(): 
                rd.fig.suptitle('_')#showme(False)
            else: 
                showme(True)
        if "holdon" in argv.keys(): 
            if "saveas" in argv.keys(): 
                gcf().savefig(argv['saveas'][0])
            else:
                showme(True)
        sys.exit()

    #-------------------------------------
    elif "savefig" in argv.keys():
        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['savefig']
        if "m96" in argv.keys():
            dms = [depthmodel_from_mod96(f) for f in argv['m96']]
        else: dms = []

        def fun(h):
            h.show(rd = None,
                plim  = argv['plim'] if 'plim' in argv.keys() else None,
                vlim  = argv['vlim'] if 'vlim' in argv.keys() else None,
                zlim  = argv['zlim'] if 'zlim' in argv.keys() else None,
                vplim = argv['vplim'] if 'vplim' in argv.keys() else None,
                vslim = argv['vslim'] if 'vslim' in argv.keys() else None,
                rhlim = argv['rhlim'] if 'rhlim' in argv.keys() else None,
                prlim = argv['prlim'] if 'prlim' in argv.keys() else None,
                nstep = int(argv['nstep'][0]) if 'nstep' in argv.keys() else 1)
            for dm in dms: 
                rd.plotmodel(dm.vp.z, dm.vp.values, dm.vs.values, dm.rh.values, color ="k", linewidth = 3, alpha = 1.)

            png = h.hlf.replace('.hlf', '.png')
            gcf().savefig(png)
            plt.close(gcf())
            #print "%s -> %s" % (h.hlf, png)

        def gen():
            for hlf in hlfs:
                h = HerrLinFile(hlf)
                yield Job(h)

        wb = waitbar('progress (figures)')
        with MapAsync(fun, gen(), Nworkers = nworkers) as ma:
            for n, _ in enumerate(ma):
                wb.refresh(n / float(len(hlfs)))
            wb.close()
        sys.exit()
    #-------------------------------------
    elif "merge" in argv.keys():
        def merge_hlfs(hlfs):
            dms = []
            zmax = -np.inf
            vsmin, vsmax = np.inf, -np.inf
            vpmin, vpmax = np.inf, -np.inf
            rhmin, rhmax = np.inf, -np.inf
            for hlf in hlfs:
                h = HerrLinFile(hlf)
                dm = depthmodel_from_mod96_txt(h.models[-1])
                dms.append(dm)
                zmax  = np.max([zmax, 1.1 * dm.vs.z[-1]])
                vsmin = np.min([vsmin, dm.vs.values.min()])
                vsmax = np.max([vsmax, dm.vs.values.max()])
                vpmin = np.min([vpmin, dm.vp.values.min()])
                vpmax = np.max([vpmax, dm.vp.values.max()])
                rhmin = np.min([rhmin, dm.rh.values.min()])
                rhmax = np.max([rhmax, dm.rh.values.max()])
            vspdf = depthpdf(z = np.linspace(0., zmax, 100), v = np.linspace(vsmin, vsmax, 100)) 
            vppdf = depthpdf(z = np.linspace(0., zmax, 100), v = np.linspace(vpmin, vpmax, 100))
            rhpdf = depthpdf(z = np.linspace(0., zmax, 100), v = np.linspace(rhmin, rhmax, 100))
            for dm in dms:
                vspdf.append(dm.vs)
                vppdf.append(dm.vp)
                rhpdf.append(dm.rh)
            zmed, vsmed = vspdf.purcentile(0.5)
            _,    vpmed = vppdf.purcentile(0.5)
            _,    rhmed = rhpdf.purcentile(0.5)

            _, vs16 = vspdf.purcentile(0.16)
            _, vs84 = vspdf.purcentile(0.84)

            _, vp16 = vppdf.purcentile(0.16)
            _, vp84 = vppdf.purcentile(0.84)

            _, rh16 = rhpdf.purcentile(0.16)
            _, rh84 = rhpdf.purcentile(0.84)


            dmmed = depthmodel(\
                    depthmodel1D(zmed, vpmed),
                    depthmodel1D(zmed, vsmed),
                    depthmodel1D(zmed, rhmed)).simplify()

            dm16 = depthmodel(\
                    depthmodel1D(zmed, vp16),
                    depthmodel1D(zmed, vs16),
                    depthmodel1D(zmed, rh16)).simplify()

            dm84 = depthmodel(\
                    depthmodel1D(zmed, vp84),
                    depthmodel1D(zmed, vs84),
                    depthmodel1D(zmed, rh84)).simplify()

            return dmmed, dm16, dm84

        #assert "hlf" in argv.keys()
        #hlfs = argv['hlf']
        hlfs = argv['merge']
        #------------------
        if len(hlfs) == 1 and not os.path.exists(hlfs[0]):
            #assume path
            hlfs = np.asarray(glob.glob(hlfs[0]))
        #------------------
        roots = []
        for hlf in hlfs:
            assert "/start" in hlf and hlf.endswith('.hlf')
            roots.append(hlf.split('/start')[0].split('/')[-1])
        roots = np.array(roots)
        for root in np.unique(roots):
            I = roots == root
            ls = hlfs[I] #list of files like pix_9_9_lon7.920235_lat48.943362_/start????.hlf to merge
            for _ in ls: 
                print _
            dmmed, dm16, dm84 = merge_hlfs(ls)
            fmedout = "%s/median.mod96" % root
            f16out = "%s/p16.mod96" % root
            f84out = "%s/p84.mod96" % root
            print "    > %s" % (fmedout, )
            dmmed.write96(fmedout, overwrite = True)
            dm16.write96(f16out, overwrite = True)
            dm84.write96(f84out, overwrite = True)
        sys.exit()
    #-------------------------------------
    elif "cube" in argv.keys():
#        assert "hlf" in argv.keys()
#        hlfs = argv['hlf']
        filesin = argv['cube']

        utmzone = argv["utm"][0].strip('N').strip('S')
        if   argv["utm"][0].endswith("N"): utmns   = "north" 
        elif argv["utm"][0].endswith("S"): utmns   = "south" 
        else: raise
        #z   = np.sort(np.array(argv['dpt']))      
        z = np.arange(*argv['dpt'])

        lons = np.asarray([readfloatafter(f, "lon") for f in filesin])
        lats = np.asarray([readfloatafter(f, "lat") for f in filesin])           
        if not len(munique(lons, lats)[0]) == len(lons):
            raise Exception('several files found with same coordinates, see --merge if --init was run with -mstart')
        lon0 = np.mean(lons)
        lat0 = np.mean(lats)
        
        LON, LAT, Z, VS, VP, RH = [], [], [], [], [], []
        for f, lon, lat in zip(filesin, lons, lats):
            if f.endswith(".mod96"):
                dm = depthmodel_from_mod96(f) #not a hlf file
            elif f.endswith('.hlf'):
                h = HerrLinFile(f)
                nlast, last = len(h.models) - 1, h.models[-1]
                dm = depthmodel_from_mod96_txt(last)
            vs = dm.vs.interp(z)
            vp = dm.vp.interp(z)
            rh = dm.rh.interp(z)
            
            LON.append(lon * np.ones_like(z))
            LAT.append(lat * np.ones_like(z))
            Z.append(z)
            VS.append(vs)
            VP.append(vp)
            RH.append(rh)
        #----------------
        LON = np.concatenate(LON)
        LAT = np.concatenate(LAT)
        Z   = np.concatenate(Z)
        VS  = np.concatenate(VS)
        VP  = np.concatenate(VP)
        RH  = np.concatenate(RH)

        #----------------     
        WGS84="+init=EPSG:4326"
        UTM="+proj=utm +zone={utmzone} +{utmns} +datum=WGS84 +units=m +no_defs".format(utmzone = utmzone, utmns=  utmns)
        p = Proj(lon0, lat0, elev0=0., wref=WGS84, lref=UTM, cartunit="km")

        X, Y = p(LON, LAT)
        cVS = cube(X, Y, Z, VS)
        cVP = cube(X, Y, Z, VP)
        cRH = cube(X, Y, Z, RH)

        #save projection information
        with open('out.proj', 'w') as fid:
            fid.write("lon0=%f\n" % lon0)
            fid.write("lat0=%f\n" % lat0)
            fid.write("elev0=%f\n" % 0.)
            fid.write('wref="%s"\n' % WGS84)
            fid.write('lref="%s"\n' % UTM)
            fid.write('cartunit="%s"\n' % "km")
        #save cube
        if     np.asarray(["median.mod96" in f for f in filesin]).all(): sfx = "med"
        elif   np.asarray(["p16.mod96"    in f for f in filesin]).all(): sfx = "p16"
        elif   np.asarray(["p84.mod96"    in f for f in filesin]).all(): sfx = "p84"
        else: sfx = ""
        with open('VS%s.xyzv' % sfx, 'w') as fid:
            for xyzv in zip(cVS.X, cVS.Y, cVS.Z, cVS.V):
                fid.write('%f %f %f %f\n' % xyzv)
        with open('VP%s.xyzv' % sfx, 'w') as fid:
            for xyzv in zip(cVP.X, cVP.Y, cVP.Z, cVP.V):
                fid.write('%f %f %f %f\n' % xyzv)
        with open('RH%s.xyzv' % sfx, 'w') as fid:
            for xyzv in zip(cRH.X, cRH.Y, cRH.Z, cRH.V):
                fid.write('%f %f %f %f\n' % xyzv)
        #----------------
        if "isec" in argv.keys(): #interactive sections
            dpt1 = float(argv['isec'][0])
            dpt2 = float(argv['isec'][1])
            dpt3 = float(argv['isec'][2])
            plt.figure()
            ax1 = gcf().add_subplot(3, 3, 1, aspect = 1.0); ax1 = gca()
            ax2 = gcf().add_subplot(3, 3, 2, aspect = 1.0, sharex = ax1, sharey = ax1); ax2 = gca()
            ax3 = gcf().add_subplot(3, 3, 3, aspect = 1.0, sharex = ax1, sharey = ax1); ax3 = gca()
            ax4 = gcf().add_subplot(3, 1, 2); ax4 = gca()
            ax5 = gcf().add_subplot(3, 1, 3); ax5 = gca()
            ax4.invert_yaxis()
            ax5.invert_yaxis()

            x = np.linspace(cVS.X.min(), cVS.X.max(), 50)
            y = np.linspace(cVS.Y.min(), cVS.Y.max(), 54)
            X, Y = np.meshgrid(x, y)
            vmin, vmax = cVS.minmax()
            cmap = tomocmap1()#plt.cm.jet
            cmap = lartceps()
#            dpt1 = z.min() + 0.1 * (z.max() - z.min())
#            dpt2 = z.min() + 0.5 * (z.max() - z.min())
#            dpt3 = z.min() + 0.9 * (z.max() - z.min())
            for z, ax in zip([dpt1, dpt2, dpt3], [ax1, ax2, ax3]):
                Z = cVS.slice(X, Y, z, method = "linear")

                if "hsmth" in argv.keys():
                    I = np.isnan(Z)
                    Z[I] = np.median(Z[~I])
                    Z, _ = gaussianblur(Z, int(argv['hsmth'][0]), int(argv['hsmth'][0]))
                    Z[I] = np.nan

                Z = np.ma.masked_where(np.isnan(Z), Z)
                ax.contourf1(X, Y, Z, vmin = vmin, vmax = vmax, cmap = cmap)
                ax.set_title('z = %.2f' % z)

            #cmap = cimseis()
            cVS.interactive_section(axsection = ax4, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "A", vmin = vmin, vmax = vmax, cmap = cmap, levels = np.linspace(vmin, vmax, 24))
            cVS.interactive_section(axsection = ax5, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "B", vmin = vmin, vmax = vmax, cmap = cmap, levels = np.linspace(vmin, vmax, 24))

            showme()
               
        sys.exit()
    #-------------------------------------
    elif "showxyzv" in argv.keys(): #HerrLin4_1 --showxyzv VSmed.xyzv -dpt 0. 3. 0.05 -isec 0.2 0.5 2.0
        c = cube_from_xyzv(argv['showxyzv'][0])


        dpt1 = float(argv['isec'][0])
        dpt2 = float(argv['isec'][1])
        dpt3 = float(argv['isec'][2])
        plt.figure()
        ax1 = gcf().add_subplot(3, 3, 1, aspect = 1.0); ax1 = gca()
        ax2 = gcf().add_subplot(3, 3, 2, aspect = 1.0, sharex = ax1, sharey = ax1); ax2 = gca()
        ax3 = gcf().add_subplot(3, 3, 3, aspect = 1.0, sharex = ax1, sharey = ax1); ax3 = gca()
        ax4 = gcf().add_subplot(3, 1, 2); ax4 = gca()
        ax5 = gcf().add_subplot(3, 1, 3); ax5 = gca()
        ax4.invert_yaxis()
        ax5.invert_yaxis()

        x = np.linspace(c.X.min(), c.X.max(), 50)
        y = np.linspace(c.Y.min(), c.Y.max(), 54)
        X, Y = np.meshgrid(x, y)
        vmin, vmax = c.minmax()
        cmap = tomocmap1()#plt.cm.jet
        cmap = lartceps()
        #cmap = toh()

        for z, ax in zip([dpt1, dpt2, dpt3], [ax1, ax2, ax3]):
            Z = c.slice(X, Y, z, method = "linear")

            if "hsmth" in argv.keys():
                I = np.isnan(Z)
                Z[I] = np.median(Z[~I])
                Z, _ = gaussianblur(Z, int(argv['hsmth'][0]), int(argv['hsmth'][0]))
                Z[I] = np.nan

            Z = np.ma.masked_where(np.isnan(Z), Z)
            ax.contourf1(X, Y, Z, vmin = vmin, vmax = vmax, cmap = cmap)
            ax.plot(c.X, c.Y, "k+", alpha = 0.1)
            ax.set_title('z = %.2f' % z)

        #cmap = cimseis()
        c.interactive_section(axsection = ax4, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "A", vmin = vmin, vmax = vmax, cmap = cmap, levels = np.linspace(vmin, vmax, 24))
        c.interactive_section(axsection = ax5, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "B", vmin = vmin, vmax = vmax, cmap = cmap, levels = np.linspace(vmin, vmax, 24))

        showme()
               
        sys.exit()



















        
