#from labex.tools.generic_nobspy import count, meshgrid3, dichotomizator
from labex.tools.generic_nobspy import meshgrid3_1
from labex.graphictools.gutils import makecolorbar, grapher, Axes3D
from labex.graphictools.pickers import pickermultiax
#from labex.graphictools.carto import mapper
#from labex.graphictools.pickers import pickermultiax
from scipy.interpolate import griddata
from labex.signaltools.smoothing import gaussianblur
from labex.signaltools.hallebach import hallebach3D
import matplotlib.pyplot as plt
import numpy as np
import copy
from labex import * #debugging

"""
regular grids (i.e. grids that can be determined by three 1D arrays with distinct lenghts) 
can be problematic in case of coordinate conversions
e.g. a regular grid in longitude, latitude, depth is no longer regular once converted into cartesian coordinates

here, we build a class that can handle 3D grids with varying coordinates in each cell of the grid
the coordinates must be expressed in a normed and orthogonal coordinate system

assume you have a regular grid given by
lon = 1d, Nlon items
lat = 1d, Nlat items
dpt = 1d, Ndpt items
V   = 3d, Nlon * Nlat * Ndpt items

1) recover the 3d LON, LAT, DPT arrays using appropriate meshgrid (depending on how V is indexed)
2) convert LON, LAT arrays into X, Y cartesian ones using for instance pyproj or labex.geotools.myproj.Proj, make sure X, Y, DPT are expressed with same units
3) initiate cube with X, Y, DPT, V
4) make slices, sections, ...
5) eventually, make backward projections to display data as a function of longitude, latitude

"""

def is3D(ax):    return hasattr(ax, 'get_proj')
def is2D(ax):    return not is3D(ax)
def check3d(ax): assert is3D(ax)
def check2d(ax): assert not is3D(ax)
def isarr(w):    return hasattr(w, "__iter__")
#--------------------------
class cube(object):
    def __init__(self, X, Y, Z, V = None):
        "indexation order does not matter"
        assert X.shape == Y.shape == Z.shape
        if V is not None: assert V.shape == Z.shape
        else: V = np.zeros_like(Z) * np.nan
        self.X = X.flat[:]
        self.Y = Y.flat[:]
        self.Z = Z.flat[:]
        self.V = V.flat[:]
    #--------------------------
    def blur(self, thickness):
        r = np.arange(len(self.X))
        for x, y in zip(*munique(self.X, self.Y)):
            I = r[(self.X == x) & (self.Y == y)]
            I = I[np.argsort(self.Z[I])]

            z, v = self.Z[I], self.V[I]
            dm = depthmodel1D(z, v)
            dm.blur(thickness)
            self.V[I] = dm.values

    #--------------------------
    def order(self, precision = 2):
        "order data so that x is the most rapidely changing value, then y and finally z"
        I = np.lexsort((np.round(self.X, precision), np.round(self.Y, precision), np.round(self.Z, precision)))
        self.X, self.Y, self.Z, self.V = [_[I] for _ in self.X, self.Y, self.Z, self.V]
    #--------------------------
    def regularize(self, x0, y0, z0, nx, ny, nz, dx, dy, dz, lambdaxmin, lambdaymin, lambdazmin, order=4., Niter=5):
        """use a hallbach algorithm to restore a regular grid using 3D lowpass butterworth filtering
           might also be used to smooth the data in 3D
        """
        x = x0 + np.arange(nx) * dx
        y = y0 + np.arange(ny) * dy
        z = z0 + np.arange(nz) * dz
        V = hallebach3D(
            xi = self.X.flat[:], 
            yi = self.Y.flat[:], 
            zi = self.Z.flat[:], 
            vi = self.V.flat[:],
            x  = x, y = y, z = z,
            fxmax = None if lambdaxmin is None else 1. / lambdaxmin, 
            fymax = None if lambdaymin is None else 1. / lambdaymin, 
            fzmax = None if lambdazmin is None else 1. / lambdazmin, 
            order=order, 
            Niter=Niter)
        X, Y, Z = meshgrid3_1(x, y, z)
        self.X = X.flat[:]
        self.Y = Y.flat[:]
        self.Z = Z.flat[:]
        self.V = V.flat[:]
    #--------------------------
    def minmax(self):
        I = np.isnan(self.V) | np.isinf(self.V)
        return minmax(self.V[~I])
    #--------------------------
    def showgrid(self, ax):
        check3d(ax)
        ax.plot(self.X, self.Y, self.Z, 'k+')
    #--------------------------
    def scatter(self, ax, **kwargs):
        check3d(ax)
        ax.scatter(self.X, self.Y, self.Z, s = 20. * np.ones_like(self.X), c = self.V, **kwargs)
    #--------------------------
    def whereami(self, x, y, z, check = True):
        """assume a normed and orthogonal coordinate system"""
        if check:
            assert not hasattr(x, "__iter__")
            assert not hasattr(y, "__iter__")
            assert not hasattr(z, "__iter__")
        d = np.sqrt((self.X - x) ** 2. + \
                    (self.Y - y) ** 2. + \
                    (self.Z - z) ** 2.)
        return np.argmin(d)
    #--------------------------
    def wherearewe(self, x, y, z):
        assert x.shape == y.shape == z.shape
        return np.asarray([self.whereami(xx, yy, zz, check = False) for xx, yy, zz in zip(x.flat[:], y.flat[:], z.flat[:])])
    #--------------------------
    def slice(self, X, Y, Z, method='linear'):
        points = np.concatenate(([self.X], [self.Y], [self.Z]), axis = 0).T
        values = self.V
        if isarr(X) and isarr(Y) and not isarr(Z):
            #horizontal slice
            assert X.shape == Y.shape
            return griddata(points, values, (X, Y, Z * np.ones_like(X)), method = method)
        elif isarr(X) and not isarr(Y) and isarr(Z):
            #west-east section
            assert X.shape == Z.shape
            return griddata(points, values, (X, Y * np.ones_like(X), Z), method = method)
        elif not isarr(X) and isarr(Y) and isarr(Z):
            #south-north section
            assert Y.shape == Z.shape
            return griddata(points, values, (X * np.ones_like(Y), Y, Z), method = method)
        else:
            raise Exception('expects 2 arrays and one float')
    #---------------------
    def slicealong(self, xpath, ypath, z, N = 100, method = "linear"):
        """
        input
            xpath, ypath   : 1d, coordinates of the path
            z              : 1d, depth array
        output
            l              : 1d, length along section
            z              : same as input
            V              : 2d, section data
            xs             : 1d, x coordinates corresponding to l
            ys             : 1d, y coordinates corresponding to l
            lpath          : 1d, length along section corresponding to input xpath and ypath
        """
        assert len(xpath) == len(ypath) >= 2

        assert len(xpath) == len(ypath)
        d = ((xpath[1:] - xpath[:-1]) ** 2. + \
             (ypath[1:] - ypath[:-1]) ** 2.) ** 0.5
        lpath = np.concatenate(([0.], d.cumsum()))
        l = np.linspace(0., d.sum(), N)
        xs = np.interp(l, xp = lpath, fp = xpath)
        ys = np.interp(l, xp = lpath, fp = ypath)
        
        Xs, Z = np.meshgrid(xs, z)
        Ys, _ = np.meshgrid(ys, z)

        #---------------------
        points = np.concatenate(([self.X], [self.Y], [self.Z]), axis = 0).T
        values = self.V
        V = griddata(points, values, (Xs, Ys, Z), method = method)
        
        #---------------------
        return l, z, V, xs, ys, lpath       
    #---------------------
    def zslice3D(self, ax, X, Y, z0, method = 'linear', *args, **kwargs):
        check3d(ax)
        V = self.slice(X, Y, z0, method = method)
        ax.contour3D(X, Y, V, zdir = "z", offset = z0, *args, **kwargs)
    #---------------------
    def xslice3D(self, ax, x0, Y, Z, method = 'linear', *args, **kwargs):
        check3d(ax)
        V = self.slice(x0, Y, Z, method = method)
        ax.contour3D(V, Y, Z, zdir = "x", offset = x0, *args, **kwargs)
    #---------------------
    def yslice3D(self, ax, X, y0, Z, method = 'linear', *args, **kwargs):
        check3d(ax)
        V = self.slice(X, y0, Z, method = method)
        ax.contour3D(X, V, Z, zdir = "y", offset = y0, *args, **kwargs)
    #---------------------
    def contour3D(self, ax, x, y, z, method = 'linear', *args, **kwargs):
        check3d(ax)
        assert len(x.shape) == 1
        assert len(y.shape) == 1
        assert len(z.shape) == 1

        for xs in x:
            YY, ZZ = np.meshgrid(y, z)
            self.xslice3D(ax, xs, YY, ZZ, *args, **kwargs)#vmin = -1., vmax = 1., levels = [0.])
        for zs in z:
            XX, YY = np.meshgrid(x, y)
            self.zslice3D(ax, XX, YY, zs, *args, **kwargs)#, vmin = -1., vmax = 1., levels = [0.])
        for ys in y:
            XX, ZZ = np.meshgrid(x, z)
            self.yslice3D(ax, XX, ys, ZZ, *args, **kwargs)#, vmin = -1., vmax = 1., levels = [0.])
        ax.set_xlim(minmax(x))
        ax.set_ylim(minmax(y))
        ax.set_zlim(minmax(z))
    #---------------------
    def interactive_section(self, axsection, axmaps, mapmarker = "ko-", move = True, method = "linear", pointnames = "A", Nz = 30, Nl = 100, *args, **kwargs):
        """let user pick the coordinates of a section accross several horizontal slices
           axmaps must be a list of axes (not mappers, otherwise picks wont appear on the plot, use ax attribute of mapper object if needed)
        """

        for axmap in axmaps:
            axmap.figure.show()
            axmap.figure.canvas.draw()

        p = pickermultiax(axmaps, marker = mapmarker, move = move)
        p.connect()
        while p.connected:
            ans = raw_input('picker is still connected')           
            if ans.lower() in ['q', 'ok', 'bye', 'leave', 'exit', 'ciao', 'disconnect', 'break']:
                #fucking mac has no MB2
                p.disconnect()
        keys = ["%s%d" % (pointnames, i) for i in xrange(len(p.xs))]
        z = np.linspace(self.Z.min(), self.Z.max(), Nz)
        l, z, V, xs, ys, lpath = self.slicealong(np.array(p.xs), np.array(p.ys), z, N = Nl, method = method)
        V = np.ma.masked_where(np.isnan(V) | np.isinf(V), V)

        axsection.contourf1(l, z, V, *args, **kwargs)
       

        for xx, yy, ll, kk in zip(p.xs, p.ys, lpath, keys):
            for axmap in axmaps:
                axmap.text(xx, yy, kk, color = "k")
            axsection.plot([ll, ll], minmax(z), 'k')
            axsection.text(ll, max(z), kk, color = 'k')

        for axmap in axmaps:
            axmap.figure.show()
            axmap.figure.canvas.draw()
        axsection.figure.canvas.draw()
    #---------------------
    def writexyzv(self, filename):
        with open(filename, 'w') as fid:
            for x, y, z, v in zip(self.X, self.Y, self.Z, self.V):
                fid.write('%f %f %f %f\n' % (x, y, z, v))

#--------------------------
class cube_from_xyzv(cube):
    def __init__(self, filename):
        A = np.loadtxt(filename)
        cube.__init__(self, A[:, 0], A[:, 1], A[:, 2], A[:, 3])

#--------------------------
if __name__ == "__main__":
    if False:
        if False:
            x = np.linspace(0., 1., 10)
            y = np.linspace(10., 11., 12)
            z = np.linspace(100., 111., 13)
            IZ, IY, IX = np.mgrid[:len(z), :len(y), :len(x)]
            x, y, z = x[IX], y[IY], z[IZ]
        elif False:
            x = np.linspace(0., 1., 10)
            y = np.linspace(10., 11., 12)
            z = np.linspace(100., 111., 13)
            x, y, z = meshgrid3(x, y, z)
        elif True:
            lon = np.linspace(7.8, 8.0, 10)
            lat = np.linspace(48.8, 49.1, 12)
            z   = np.linspace(-3., 0., 13)
            mp  = MyProj(cartunit = "km") #conversion induce a small rotation
            lon, lat, z = meshgrid3(lon, lat, z)
            x, y = mp(lon, lat)
            
            

        c = cube(x, y, z)
        c.order(precision = 6)
        c.V = np.arange(len(x.flat[:]))

        plt.figure()
#        ax = Axes3D(gcf())
#        gcf().show()
        for  x, y, z, v in zip(c.X, c.Y, c.Z, c.V):
            print "%10f %10f %10f %10f" % (x, y, z, v)
#            ax.plot([x], [y], [z], 'k+')
#        ax.figure.canvas.draw()


        c.scatter(Axes3D(gcf()))
        showme()

    elif False:
        x = np.random.rand(12500) * 4. - 2.
        y = np.random.rand(12500) * 4. - 2.
        z = np.random.rand(12500) * 4. - 2.

        D = np.sqrt((x - 1.) ** 2. + (y + 2.) ** 2. + z ** 2.)
        v = np.sin(3.3 * D)

        D = np.sqrt((x + 1.) ** 2. + (y + 0.) ** 2. + (z - 3.) ** 2.)
        v += np.sin(3.3 * D)

        c = cube(x, y, z, v)  
        vmin, vmax = c.minmax()
        cmap = plt.cm.jet

        if False:
            c.scatter(Axes3D(gcf()))

        if False:
            plt.figure()
            x = np.linspace(-2., 2., 12)
            y = np.linspace(-2., 2., 12)
            z = np.linspace(-2., 2., 12)
            c.contour3D(Axes3D(gcf()), x, y, z, vmin = vmin, vmax = vmax, cmap = cmap, levels = [-0.8, +0.8], method = "linear")
            gcf().show()

        if True:
            plt.figure()
            ax1 = gcf().add_subplot(3, 3, 1); ax1 = gca()
            ax2 = gcf().add_subplot(3, 3, 2, sharex = ax1, sharey = ax1); ax2 = gca()
            ax3 = gcf().add_subplot(3, 3, 3, sharex = ax1, sharey = ax1); ax3 = gca()
            ax4 = gcf().add_subplot(3, 1, 2); ax4 = gca()
            ax5 = gcf().add_subplot(3, 1, 3); ax5 = gca()

            x = np.linspace(-2., 2., 54)
            y = np.linspace(-2., 2., 65)
            X, Y = np.meshgrid(x, y)



            for z, ax in zip([-1.5, 0., 1.0], [ax1, ax2, ax3]):
                Z = c.slice(X, Y, z, method = "linear")
                Z = np.ma.masked_where(np.isnan(Z), Z)
                ax.contourf1(X, Y, Z, vmin = vmin, vmax = vmax, cmap = cmap)
                ax.set_title('z = %.2f' % z)


            c.interactive_section(axsection = ax4, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "A", vmin = vmin, vmax = vmax, cmap = cmap)
            c.interactive_section(axsection = ax5, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "B", vmin = vmin, vmax = vmax, cmap = cmap)

        showme()

    elif True:
        x = np.random.rand(12500) * 4. - 2.
        y = np.random.rand(12500) * 4. - 2.
        z = np.random.rand(12500) * 4. - 2.

        D = np.sqrt((x - 1.) ** 2. + (y + 2.) ** 2. + z ** 2.)
        v = np.sin(3.3 * D)

        D = np.sqrt((x + 1.) ** 2. + (y + 0.) ** 2. + (z - 3.) ** 2.)
        v += np.sin(3.3 * D)

        v += np.random.randn(len(v))

        c = cube(x, y, z, v)  

        c.regularize(x0 = -2., y0 = -2., z0 = -2., 
            nx = 12, ny = 14, nz = 16, 
            dx = 4. / 12., dy = 4. / 14., dz = 4. / 16., 
            lambdaxmin = 1., 
            lambdaymin = 1., 
            lambdazmin = 1., 
            order=4., Niter=5)

        vmin, vmax = c.minmax()
        cmap = plt.cm.jet


        if True:
            plt.figure()
            ax1 = gcf().add_subplot(3, 3, 1); ax1 = gca()
            ax2 = gcf().add_subplot(3, 3, 2, sharex = ax1, sharey = ax1); ax2 = gca()
            ax3 = gcf().add_subplot(3, 3, 3, sharex = ax1, sharey = ax1); ax3 = gca()
            ax4 = gcf().add_subplot(3, 1, 2); ax4 = gca()
            ax5 = gcf().add_subplot(3, 1, 3); ax5 = gca()

            x = np.linspace(-2., 2., 54)
            y = np.linspace(-2., 2., 65)
            X, Y = np.meshgrid(x, y)



            for z, ax in zip([-1.5, 0., 1.0], [ax1, ax2, ax3]):
                Z = c.slice(X, Y, z, method = "linear")
    
#                print "smoothing"
#                I = np.isnan(Z)
#                Z[I] = np.median(Z[~I])
#                Z, _ = gaussianblur(Z, 3, 3)
#                Z[I] = np.nan

                Z = np.ma.masked_where(np.isnan(Z), Z)
                ax.contourf1(X, Y, Z, vmin = vmin, vmax = vmax, cmap = cmap)
                ax.set_title('z = %.2f' % z)


            c.interactive_section(axsection = ax4, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "A", vmin = vmin, vmax = vmax, cmap = cmap)
            c.interactive_section(axsection = ax5, axmaps = [ax1, ax2, ax3], mapmarker = "ko-", move = True, method = "linear", pointnames = "B", vmin = vmin, vmax = vmax, cmap = cmap)

        showme()

            
