#!/usr/bin/env python

help = '''m96
--show            list of mod96 files to display (same plot)
--disp            name of mod96 file to use as input
    -RU0          rayleigh, group, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    -RU1          rayleigh, group, mode 1 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    -RC0          rayleigh, phase, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale          
    -LC0          love,     phase, mode 0 : expects 4 frequency arguments : fstart, fend, nfreq, fscale
    ...
    -save         name of surf96file to write
--split 
    -thck         thickness of the sublayers in km
    -sfx          suffix to add before file extension (default split)
    -o            ignore suffix and overwrite input file
--addlayer 
    -thck         thickness of the sublayers in km
    -sfx          suffix to add before file extension (default split)
    -o            ignore suffix and overwrite input file    
-inline           replace showme by plt.show (e.g. for jupyter)
'''

# ---------------------------------------------
if __name__ == "__main__":
    import sys, os
    if len(sys.argv) == 1:
        print help
        sys.exit()

    from srfpython.standalone.display import gcf, gca, showme, pause, plt, logtick
    from srfpython.utils import readargv
    from srfpython.Herrmann.Herrmann import dispersion_1
    from srfpython.depthdisp.dispcurves import surf96reader, Claw, freqspace
    from srfpython.depthdisp.depthmodels import depthmodel_from_mod96



    import numpy as np
    argv = readargv()

    # -----------------------------------
    if "inline" in argv.keys():
        showme = plt.show

    # -----------------------------------
    if "help" in argv.keys() or "h" in argv.keys():
        print help
        sys.exit()
    
    # -----------------------------------
    elif "split" in argv.keys():
        thickness = argv['thck'][0]
        suffix    = argv['sfx'][0] if "sfx" in argv.keys() else "split"

        for f in argv['split']:
            dm = depthmodel_from_mod96(f)
            dm.split(thickness)
            if "o" in argv.keys():
                fout = f
            else:
                fout = f.split('.')[:-1] + [suffix,  f.split('.')[-1]]

                fout = ".".join(fout)
                assert not os.path.exists(fout)
            print (f, ">", fout)
            dm.write96(fout)
        sys.exit()

    # -----------------------------------
    elif "addlayer" in argv.keys():
        thickness = argv['thck'][0]
        suffix    = argv['sfx'][0] if "sfx" in argv.keys() else "add"
        for f in argv['addlayer']:
            dm = depthmodel_from_mod96(f)
            dm.add_layer(thickness)
            if "o" in argv.keys():
                fout = f
            else:
                fout = f.split('.')[:-1] + [suffix,  f.split('.')[-1]]
                fout = ".".join(fout)
                assert not os.path.exists(fout)
            print (f, ">", fout)
            dm.write96(fout)
        sys.exit()

    # -----------------------------------
    elif "show" in argv.keys():
        axvp = gcf().add_subplot(1, 4, 1)
        axvs = gcf().add_subplot(1, 4, 2, sharey = axvp)
        axpr = gcf().add_subplot(1, 4, 3, sharey = axvp)
        axrh = gcf().add_subplot(1, 4, 4, sharey = axvp)
        for f in argv['show']:
            dm = depthmodel_from_mod96(f)
            print f

            dm.vp.show(axvp, ".-")
            dm.vs.show(axvs, ".-")
            dm.pr().show(axpr, ".-")
            dm.rh.show(axrh, ".-")

        if "gardner74" in argv.keys():
            rh = dm.rh.copy()
            rh.values = 1.74 * dm.vp.values ** 0.25
            rh.show(axrh)

        showme()
        sys.exit()

    # -----------------------------------
    elif "disp" in argv.keys():


        for m in argv['disp']:
            dm = depthmodel_from_mod96(m)
            ztop = dm.vp.ztop()
            vp   = dm.vp.values
            vs   = dm.vs.values
            rh   = dm.rh.values

            Waves, Types, Modes, Freqs = [], [], [], []
            for k in argv.keys():
                if k[0].upper() in "RL" and k[1].upper() in "UC" and k[2] in "0123456789":
                    fstart, fend, nfreq, fspace = argv[k]
                    freq = freqspace(float(fstart), float(fend), int(nfreq), fspace)
                    Waves.append(k[0])
                    Types.append(k[1])
                    Modes.append(int(k[2:]))
                    Freqs.append(freq)

            if "save" in argv.keys():
                sfx = ""
                s96out = ".".join(m.split('/')[-1].split('.')[:-1]) + sfx + ".surf96"
                while os.path.exists(s96out):
                    s96out = s96out.split(sfx + '.surf96')[0]
                    if sfx == "": sfx = "_1"
                    else: sfx = "_%d" % (int(sfx.strip("_")) + 1)
                    s96out = s96out + sfx + ".surf96"
                print "%s => %s" % (m, s96out)

                assert s96out.endswith('.surf96') or s96out.endswith('.s96')
                with open(s96out, 'w') as fid:
                    for w, t, m, F, V in dispersion_1(ztop, vp, vs, rh, Waves, Types, Modes, Freqs, h=0.0005, dcl=0.00005, dcr=0.00005):
                        for FF, VV in zip(F, V):
                            fid.write('SURF96 %s %s T %d %f %f 0.1\n' % (w, t, m, 1. / FF, VV))
            else:
                axvp = gcf().add_subplot(1, 5, 1)
                axvs = gcf().add_subplot(1, 5, 2, sharey=axvp)
                axpr = gcf().add_subplot(1, 5, 3, sharey=axvp)
                axrh = gcf().add_subplot(1, 5, 4, sharey=axvp)
                axdsp = gcf().add_subplot(2, 5, 5)
                dm.vp.show(axvp)
                dm.vs.show(axvs)
                dm.pr().show(axpr)
                dm.rh.show(axrh)

                for w, t, m, F, V in dispersion_1(ztop, vp, vs, rh, Waves, Types, Modes, Freqs):
                    axdsp.loglog(1. / F, V, label="%s%s%d" % (w, t, m))

                logtick(axdsp, "xy")
                plt.legend()
                showme()
        sys.exit()



















        
