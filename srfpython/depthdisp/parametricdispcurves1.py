import numpy as np
from matplotlib import ticker


"""
create a arctan-shaped phase dispersion curve in the log(p) - log(c) domain with 4 parameters
"""

def btan(x):
    return np.arctan(x) / np.pi + 1./2.

class PhaseDispersionLaw(object):
    def __init__(self, cmin=2., cmax=3.0, sharpness=1.0, fmid=1.0):
        """
        :param cmin: lower phase velocity km/s
        :param cmax: upper phase velocity km/s
        :param sharpness: sharpness of the dispersion curve at fmid
        :param fmid: inflection frequency (Hz)
        """
        self.cmin = cmin
        self.cmax = cmax
        self.sharpness = sharpness
        self.fmid = fmid
        self.eps = 1.e-20
        
    def __call__(self, freq):
        """
        evaluate the phase dispersion law at freq (array, Hz)
        """
        #return self.cmin * (self.cmax / self.cmin) ** (.5 + np.arctan(self.sharpness * np.log(self.fmid / (np.abs(freq) + self.eps))) / np.pi)
        return np.exp(np.log(self.cmin) + np.log(self.cmax / self.cmin) * (btan(self.sharpness * (np.log(self.fmid / (np.abs(freq) + self.eps))))))

    def expr1(self):
        fmt = r"$c(f) = exp \left( ln%s +  ln \frac{%s}{%s}  \left( { \frac{1}{2} + \frac{1}{\pi} atan\left(%s ln \frac{%s}{|f|}\right) } \right) \right)$"
        return fmt % ("{}".format(self.cmin),
                      "{}".format(self.cmax),
                      "{}".format(self.cmin),
                      "{}".format(self.sharpness) if self.sharpness != 1.0 else "",
                      "{}".format(self.fmid))

    def expr2(self):
        fmt = r"$c(f) =  %s \left( \frac{%s}{%s} \right) ^{ \frac{1}{2} + \frac{1}{\pi} atan(%s ln \frac{%s}{|f|}) } $"
        return fmt % ("{}".format(self.cmin),
                      "{}".format(self.cmax),
                      "{}".format(self.cmin),
                      "{}".format(self.sharpness) if self.sharpness != 1.0 else "",
                      "{}".format(self.fmid))

    def show(self, ax, freqmin, freqmax, **kwargs):
        f = np.logspace(np.log10(freqmin), np.log10(freqmax), 100)
        ax.loglog(1. / f, self(f), **kwargs)

        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_locator(ticker.LogLocator(base = 10., subs=[1., 2., 3., 5.]))
            axis.set_minor_locator(ticker.NullLocator())
            axis.set_major_formatter(ticker.StrMethodFormatter('{x}'))
        ax.set_xlabel('period (s)')
        ax.set_ylabel('velocity (km/s)')

    def to_ulaw(self):
        return GroupDispersionLaw(self)

class GroupDispersionLaw(PhaseDispersionLaw):
    def __init__(self, phasedispersionlaw):
        self.pdl = phasedispersionlaw
        self.df = 0.01

    def expr1(self):
        fmt = r"$ u(f) = c(f) / (1. - \frac{ln\,c(f+ \delta f) - ln\,c(f)}{ln(f+\delta f)-ln(f)})$"
        return fmt

    def __call__(self, freq):
        return self.pdl(freq) / (1. - (np.log(self.pdl(freq+self.df) / self.pdl(freq)) / np.log((freq+self.df)/freq)))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    c = PhaseDispersionLaw(cmin=2.25, cmax=3.52, sharpness=1.1, fmid=0.25)
    u = GroupDispersionLaw(c)

    c.show(plt.gca(), 0.05, 1.0, label=c.expr1())
    f = np.logspace(np.log10(0.05), 0., 100, 'plog')
    plt.gca().plot(1./f, c(f), label=c.expr2())
    plt.show()
