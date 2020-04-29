#!/usr/bin/env python

help = '''s96
--show            list of surf96files to display
    -freq         plot in frequency domain if specified
--resamp          list of surf96files to resample
    -fspace       new frequency array in Hz, fstart, fend, nfreq, fscale
    -sfx          file suffix to append, use "" to overwrite input files
-inline             replace showme by plt.show (e.g. for jupyter)
#surf96 format 
SURF96 {wave} {type} {flag} {mode} {period(s)} {value(km/s)} {dvalue(km/s)}
'''

# ---------------------------------------------
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        print help
        sys.exit()

    from srfpython.standalone.display import gcf, gca, showme, pause, plt, logtick
    from srfpython.depthdisp.dispcurves import surf96reader, Claw, freqspace
    from srfpython.utils import readargv
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
    elif "show" in argv.keys():
        pmin, pmax = np.inf, -np.inf
        for f in argv['show']:
            s = surf96reader(f)
            print f
            for law in s.get_all():
                print "    %s" % str(law)
                law.show(gca(), period= not "freq" in argv.keys(), label="%s%s%d" % (law.wave, law.type, law.mode))#alpha = 0.5, color = "r" if law.mode else "k")
                pmin = min([pmin, 1. / law.freq.max()])
                pmax = max([pmax, 1. / law.freq.min()])

        gca().set_xscale('log')
        gca().set_yscale('log')
        logtick(gca(), "xy")
        if not "freq" in argv.keys(): gca().set_xlim(pmin, pmax)
        else:                         gca().set_xlim(1. / pmax, 1. / pmin)
        gca().grid(True)
        gca().set_xlabel("period (s)")
        gca().set_ylabel("velocity (km/s)")
        # plt.legend()
        showme()
        sys.exit()

    # -----------------------------------
    elif "resamp" in argv.keys():
        def tostr(law, newf):
            fmt = "SURF96 {wave} {type} {flag} {mode} {period} {value} {dvalue}"
            stdlaw = Claw(freq = law.freq, value = law.dvalue, extrapolationmode = 0)
            newv = law(newf)
            news = stdlaw(newf)
            I = ~np.isnan(newv)            
            if I.any():
                s = "\n".join([fmt.format(\
                            wave = law.wave, 
                            type = law.type, 
                            flag = law.flag, 
                            mode = law.mode, 
                            period = 1./ff, 
                            value  = vv, 
                            dvalue = ss) for ff, vv, ss in zip(newf[I], newv[I], news[I])])
                return s
            return ""

        sfx  = argv["sfx"][0] if "sfx" in argv.keys() else "resamp"        
        print argv['fspace']
        newf = freqspace(float(argv['fspace'][0]), float(argv['fspace'][1]), int(argv['fspace'][2]), argv['fspace'][3]) #np.sort(np.unique(np.abs(argv['newf'])))
        print newf
        for f in argv['resamp']:
            s = surf96reader(f)
            fout = ".".join(f.split('/')[-1].split('.')[:-1] + [sfx] + [f.split('/')[-1].split('.')[-1]])
            print fout
            out = ""
            for law in s.get_all():
                out += "\n" + tostr(law, newf)
            out = out.strip('\n').strip()
            with open(fout, 'w') as fid:
                fid.write(out)


            


