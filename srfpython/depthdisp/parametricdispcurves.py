import numpy as np
from matplotlib import ticker


"""
create a arctan-shaped phase dispersion curve in the log(p) - log(c) domain with 4 parameters
"""

def btan(x):
    return np.arctan(x) / np.pi + 1./2.

class PhaseDispersionLaw(object):
    def __init__(self, cmin=2., cmax=3.0, power=1.0, pmid=1.0):
        """
        :param cmin: lower phase velocity km/s
        :param cmax: upper phase velocity km/s
        :param power: sharpness of the dispersion curve at pmid
        :param pmid: inflection period (s)
        """
        self.cmin = cmin
        self.cmax = cmax
        self.power = power
        self.pmid = pmid

        self.offset = np.log(cmin)
        self.gain = np.log(cmax / cmin)
        self.power = power
        self.inflection=np.log(pmid)
        self.eps = 1.e-20
        
    def __call__(self, freq):
        """
        evaluate the phase dispersion law at freq (array, Hz)
        """
        logp = np.log(1. / (np.abs(freq) + self.eps))
        logc = self.offset + self.gain * btan(self.power * (logp - self.inflection))
        return np.exp(logc)

    def tester(self, freq):
        fmid = 1. / self.pmid
        return self.cmin * (self.cmax / self.cmin) ** \
               (.5 + np.arctan(self.power * np.log(fmid / (np.abs(freq) + self.eps))) / np.pi)

    def show(self, ax, freqmin, freqmax, **kwargs):
        f = np.logspace(np.log10(freqmin), np.log10(freqmax), 100)
        ax.loglog(1. / f, self(f), **kwargs)

        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_locator(ticker.LogLocator(base = 10., subs=[1., 2., 3., 5.]))
            axis.set_minor_locator(ticker.NullLocator())
            axis.set_major_formatter(ticker.StrMethodFormatter('{x}'))
        ax.set_xlabel('period (s)')
        ax.set_ylabel('velocity (km/s)')


class GroupDispersionLaw(PhaseDispersionLaw):
    def __init__(self, phasedispersionlaw):
        self.pdl = phasedispersionlaw
        self.df = 0.01
        
    def __call__(self, freq):
        return self.pdl(freq) / (1. - (np.log(self.pdl(freq+self.df) / self.pdl(freq)) / np.log((freq+self.df)/freq)))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    c = PhaseDispersionLaw(cmin=2.25, cmax=3.52, power=1.0, pmid=4.0)
    u = GroupDispersionLaw(c)

    f = np.logspace(-1., 1., 100, 'plog')
    plt.gca().plot(1./f, c(f))
    plt.gca().plot(1./f, c.tester(f))
    plt.show()