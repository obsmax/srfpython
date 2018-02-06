import os
import numpy as np
from tetedenoeud.utils.display import value2color
from tetedenoeud.utils.cmaps import tej
from srfpython.depthdisp.depthmodels import depthmodel_from_mod96, \
    depthmodel_from_mod96string, depthmodel, depthmodel1D
from srfpython.depthdisp.dispcurves import surf96reader, \
    surf96reader_from_surf96string
from srfpython.depthdisp.depthdispdisplay import DepthDispDisplay, plt, showme



# ---------------------------------------------------
def m962str(m96):
    dm = depthmodel_from_mod96(m96)
    return dm.__str__()


# ---------------------------------------------------
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

    # -----------------------------------------------
    def __init__(self, hlf):
        assert isinstance(hlf, str)
        assert hlf.endswith('.hlf')
        self.hlf = hlf
        self.target = None  # a tuple like (targetfilename, s96str)
        self.init = None  # name of initfile
        self.models = []  # list of models (m96str), first = init
        self.forwards = []  # list of dispersion curves (s96str)
        self.chi2reds = []  # list of floats
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
                            continue  # next line already read in
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
                            continue  # next line already read in
                        elif ">forward" in l:
                            nmodel = int(l.split('>forward.')[-1].split()[0])
                            assert len(self.forwards) == nmodel
                            L = []
                            while True:
                                l = fid.readline()
                                if l == "" or l.startswith(">"): break
                                L.append(l)
                            self.forwards.append("".join(L).strip('\n'))
                            continue  # next line already read in
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
    # -----------------------------------------------
    def __str__(self):
        s = ""
        if self.target is None:
            return s
        else:
            targetfile, L = self.target
            s += ">target %s\n%s\n" % (targetfile, L)
            if not len(self.models):
                return s
            else:
                for nmodel in xrange(len(self.models)):
                    s += ">model.%03d %s\n%s\n" % (nmodel, "" if nmodel else self.init, self.models[nmodel])
                    if len(self.forwards) > nmodel:
                        s += ">forward.%03d\n%s\n" % (nmodel, self.forwards[nmodel])
                        s += ">chi2red.%03d\n%f\n" % (nmodel, self.chi2reds[nmodel])
        return s  # .strip('\n').strip()

    # -----------------------------------------------
    def set_target(self, filename, data):
        if self.target is not None:
            raise Exception('%s has already a target' % self.hlf)

        self.target = (filename, data)
        with open(self.hlf, 'a') as fout:
            fout.write('>target %s\n' % filename)
            fout.write('%s\n' % data)

    # -----------------------------------------------
    def set_init(self, m96):
        if self.target is None:
            raise Exception('%s has no target, please set it first' % self.hlf)
        if not len(self.models) == 0:
            raise Exception('%s has already an initial model' % self.hlf)

        if os.path.exists(m96):  # m96 is a file name
            L = m962str(m96)
            self.init = m96
        else:  # m96 is a file content
            L = m96
            self.init = "_"

        self.models.append(L)

        with open(self.hlf, 'a') as fout:
            fout.write('>model.000 %s\n' % self.init)
            fout.write('%s\n' % L)

    # ---------------------------------------------------
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
        assert not hasattr(VP, "__iter__")
        assert not hasattr(VS, "__iter__")
        assert not hasattr(RH, "__iter__")

        dm = depthmodel_from_mod96string(self.models[0])
        ztop, vp, vs, rh = dm.vp.z, dm.vp.values, dm.vs.values, dm.rh.values
        if ZTOP > ztop[-1]:
            ztop = np.concatenate((ztop, [ZTOP]))
            vp = np.concatenate((vp, [VP]))
            vs = np.concatenate((vs, [VS]))
            rh = np.concatenate((rh, [RH]))
        elif ZTOP == ztop[-1]:
            vp[-1] = VP
            vs[-1] = VS
            rh[-1] = RH
        else:
            i = np.searchsorted(ztop, ZTOP)
            ztop = np.concatenate((ztop[:i], [ZTOP]))
            vp = np.concatenate((vp[:i], [VP]))
            vs = np.concatenate((vs[:i], [VS]))
            rh = np.concatenate((rh[:i], [RH]))
        dm = depthmodel( \
            depthmodel1D(ztop, vp),
            depthmodel1D(ztop, vs),
            depthmodel1D(ztop, rh))
        self.models[0] = str(dm)
        self.init += "*"

        with open(self.hlf, 'w') as fout:
            fout.write(str(self))  # .rstrip().rstrip('\n'))

    # -----------------------------------------------
    def addmodel(self, ZTOP, VP, VS, RH):
        assert len(self.chi2reds) == len(self.forwards) == len(self.models)
        dm = depthmodel( \
            depthmodel1D(ZTOP, VP),
            depthmodel1D(ZTOP, VS),
            depthmodel1D(ZTOP, RH))
        self.models.append(dm.__str__().strip('\n').strip())

    # -----------------------------------------------
    def addforward(self, chi2red, waves, types, modes, freqs, values, dvalues=None):
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

    # -----------------------------------------------
    def show(self, rd=None, plim=None, vlim=None, zlim=None, vplim=None, vslim=None, rhlim=None, prlim=None, nstep=1,
             cmap=tej()):
        target = surf96reader_from_surf96string(self.target[1])

        if rd is None:
            rd = DepthDispDisplay()

        rd.plotdisp(color="k", linewidth=3, alpha=1., *target.wtmfvd())

        N = len(self.forwards)
        clrs = [value2color(n / float(N), cmap=cmap) for n in np.arange(N)]

        def subshow(i, linewidth=1, alpha=1, linestyle="-"):
            model = depthmodel_from_mod96string(self.models[i])
            model = (model.vp.z, model.vp.values, model.vs.values, model.rh.values)
            data = surf96reader_from_surf96string(self.forwards[i])
            data = data.wtmfv()

            rd.plotmodel(color=clrs[i], linewidth=linewidth, alpha=alpha, linestyle=linestyle, *model)
            rd.plotdisp(color=clrs[i], linewidth=linewidth, alpha=alpha, linestyle=linestyle, *data)
            rd.axconv.plot(np.arange(N)[i], self.chi2reds[i], 'o', color=clrs[i])

        if nstep == -1:  # last model only
            subshow(-1, linewidth=1, alpha=0.2)
        elif nstep == 0:  # first model only
            subshow(0, linewidth=1, alpha=0.2)
        elif nstep == -2:  # first and last models only
            subshow(0, linewidth=1, alpha=0.2)
            subshow(-1, linewidth=1, alpha=0.2)
        else:
            subshow(0, linewidth=1, alpha=1)
            # display intermediate steps
            for i in np.arange(0, N, nstep)[1:-1]:
                subshow(i, linewidth=1, alpha=0.4)
            subshow(-1, linewidth=3, alpha=1)

        if plim is not None: rd.set_plim(plim)
        if vlim is not None: rd.set_vlim(vlim)
        if zlim is not None: rd.set_zlim(zlim)

        if vplim is not None: rd.axvp.set_xlim(vplim)
        if vslim is not None: rd.axvs.set_xlim(vslim)
        if prlim is not None: rd.axpr.set_xlim(prlim)
        if rhlim is not None: rd.axrh.set_xlim(rhlim)

        rd.grid()
        rd.tick()
        # rd.fig.suptitle(self.target[0])
        rd.fig.suptitle(self.hlf)