import numpy as np
import matplotlib.pyplot as plt
from srfpython import depthmodel_from_arrays, Curve
from srfpython.Herrmann.Herrmann import HerrmannCaller
from srfpython.sensitivitykernels.sker17 import sker17_2
from srfpython.standalone.cmaps import tomocmap1
from srfpython.standalone.display import logtick, textonly

"""
experimental
"""

# ============
# build depth model : !!UNITS!! (km, km/s, km/s, g/cm3)
dm = depthmodel_from_arrays(
    z   =np.array([0.,      0.5,    1.5,    5.0,     10.,        ]) * 0.001,
    vp  =np.array([    3800.,  1800.,   1500.,  2100.,   3200.,  ]) * 0.001,
    vs  =np.array([    2030.,   600.,    500.,   800.,   1500.,  ]) * 0.001,
    rh  =np.array([    2400.,  2000.,   2000.,  2000.,   2300.,  ]) * 0.001,
    )

# add a layer in the top half space before depth resampling
dm.add_layer(thickness=10. * 0.001)  # in km

# split the continuous layers into smaller layers
dm.split(thickness=0.25 * 0.001)  # in km

if False:
    # display (qc)
    dm.show(plt.gca(), ".-", units='m/s/kg.m-3',)
    plt.show()

# define the expected dispersion curves
freqs = np.linspace(5., 300., 70) # Hz
curves = [
    Curve(wave="R", type="C", mode=0, freqs=freqs), 
    Curve(wave="R", type="C", mode=1, freqs=freqs),
    Curve(wave="R", type="C", mode=2, freqs=freqs),
    Curve(wave="R", type="C", mode=3, freqs=freqs),
    ]

# forward problem
hc = HerrmannCaller(curves, ddc=0.005, h=0.005)
curves = hc.call_dm(depthmodel=dm)  # overwrite with output disp curves

if False:   
     # display (qc)
    for curve_out in curves:
        plt.plot(
            curve_out.freqs, 
            curve_out.values * 1000., 
            label="%s%s%d" % (curve_out.wave, curve_out.type, curve_out.mode))

    plt.legend()
    plt.show()

# sensitivity kernels
norm=True
sker_gen = sker17_2(
    depthmodel=dm, 
    curves=curves,
    dlogz=0.001, dlogvp=0.01, dlogvs=0.01, dlogrh=0.01, norm=norm,
    h=0.005, ddc=0.005)  

# ===================    
for wave, typ, mode, freqs, dlogvadlogz, dlogvadlogvp, dlogvadlogvs, dlogvadlogrh, dlogvadlogpr in sker_gen:

    # sensitivity kernel
    z_edges = np.hstack((dm.vp.z, [1.1 * dm.vp.z[-1]])) * 1000.  # m
    z_mid = np.hstack((0.5 * (dm.vp.z[1:] + dm.vp.z[:-1]), [1.05 * dm.vp.z[-1]])) * 1000.0  # m
    f_edges = np.hstack((freqs[0] * 0.95, np.sqrt(freqs[1:] * freqs[:-1]), freqs[-1] * 1.05))

    if not norm:
        # mask half space because it integrates the sensitivity over very thick layer => overestimated sensitivity
        for _ in dlogvadlogz, dlogvadlogvp, dlogvadlogvs, dlogvadlogrh, dlogvadlogpr:
            _[-1, :] = np.nan
        dlogvadlogz, dlogvadlogvp, dlogvadlogvs, dlogvadlogrh, dlogvadlogpr = \
            [np.ma.masked_where(np.isnan(_), _) \
             for _ in [dlogvadlogz, dlogvadlogvp, dlogvadlogvs, dlogvadlogrh, dlogvadlogpr]]

    cmap = tomocmap1(w=0.01, W=0.2)  # cccfcmap3() #plt.cm.RdBu
    for M, p, q in zip([dlogvadlogz, dlogvadlogvp, dlogvadlogvs, dlogvadlogrh, dlogvadlogpr],
                       ["ln Z^{top}_i", r"ln Vp_i", "ln Vs_i", r"ln \rho _i", "ln (Vp/Vs)_i"],
                       ["lnZtop", "lnVp", "lnVs", "lnrho", "lnVpaVs"]):
                           
        # ==== display
        fig1 = plt.figure(figsize=(8,8))
        fig1.subplots_adjust(wspace=0.3)
        gs = fig1.add_gridspec(3, 3)
        
        # depth model
        ax1 = fig1.add_subplot(gs[1:, 0])
        dm.show(ax1, ".-", units='m/s/kg.m-3',)
        ax1.grid(True, linestyle=":", color="k")
        plt.legend()

        # disp curve
        ax2 = fig1.add_subplot(gs[0, 1:], ylabel='velocity (m/s)')
        for curve in curves:           
            ax2.plot(
                curve.freqs, 
                curve.values * 1000., 
                '+-', 
                label=f"{curve.wave}{curve.type}{curve.mode}")
                
        ax2.grid(True, which="major")
        ax2.grid(True, which="minor")
        if False:
            ax2.set_xscale('log')
            ax2.set_yscale('log')            
            logtick(ax2, "xy")
        plt.legend()
                               
        ax3 = fig1.add_subplot(gs[1:, 1:], sharex=ax2, sharey=ax1)
        ax3.set_title('%s%s%d' % (wave, typ, mode))

        vmax = abs(M).max()
        vmin = -vmax

        coll = plt.pcolormesh(f_edges, z_edges, M,
                              vmin=vmin, vmax=vmax, cmap=cmap)

        if False:
            if M.max() - M.min():
                levels = np.logspace(-1., 2, 10)
                levels = np.hstack((-levels[::-1], [0], levels))
                plt.contour(freqs, z_mid, M,
                            levels=levels,
                            colors="k")

        cax = fig1.add_axes((.91, .2, .01, .2))
        plt.colorbar(coll, cax=cax)

        ax3.set_xlabel('Frequency (Hz)')

        ax3.set_xlim(min(f_edges), max(f_edges))
        ax3.set_ylim(min(z_edges), max(z_edges))
        if False:
            ax3.set_xscale('log')
            logtick(ax3, "x")
            
        if norm:
            textonly(ax3, txt=r'$ \frac{H}{H_i} \, \frac{d ln%s_j}{d %s} $' % (typ, p), loc=3, fontsize=16)
        else:
            textonly(ax3, txt=r'$ \frac{d ln%s_j}{d %s} $' % (typ, p), loc=3, fontsize=16)

        if not ax3.yaxis_inverted():
            ax3.invert_yaxis()


    # ------
    plt.show()

# ===================
