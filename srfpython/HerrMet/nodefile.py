import os
import numpy as np
from srfpython.HerrMet.files import ROOTNAME, HERRMETPARAMFILE, HERRMETTARGETFILE, HERRMETEXTRACTPDFMODELFILE
from srfpython.standalone.asciifile import AsciiFile_fromstring

"""
a node file is an ascii file with 3 columns and one line header 

#met rootdir = "."
#met extract_mode = "best"
#met extract_limit = 1000
#met extract_llkmin = 0
#met extract_step = 1
#fld longitude latitude node
#unt deg       deg      -
#fmt %f        %f       %s
12.132131 45.355425434 node0001
22.132131 35.355425434 node0002
32.132131 25.355425434 node0003
...
"""


class NodeFileString(object):
    def __init__(self, string, rootdir="."):
        """
        :param string:
        """
        a = AsciiFile_fromstring(string)

        self.lons = np.asarray(a.data['longitude'], float)
        self.lats = np.asarray(a.data['latitude'], float)
        self.nodes = np.asarray(a.data['node'], str)

        # the file rootdir is relative to the file location, stack
        self.rootdir = os.path.join(rootdir, a.metadata['rootdir'])
        self.extract_mode = a.metadata['extract_mode']
        self.extract_limit = a.metadata['extract_limit']
        self.extract_llkmin = a.metadata['extract_llkmin']
        self.extract_step = a.metadata['extract_step']

        if not os.path.isdir(self.rootdir):
            raise IOError()

        if not len(self.nodes) == len(np.unique(self.nodes)):
            raise ValueError('node names are not unique')

        n = len(self.lons)
        self.paramfiles = np.zeros(n, object)  # object for now
        self.targetfiles = np.zeros(n, object)  # object for now
        self.medianfiles = np.zeros(n, object)  # object for now
        self.p16files = np.zeros(n, object)  # object for now
        self.p84files = np.zeros(n, object)  # object for now

        for i, (lon, lat, node) in enumerate(zip(self.lons, self.lats, self.nodes)):

            rootname = os.path.join(self.rootdir, ROOTNAME.format(node=node))
            paramfile = HERRMETPARAMFILE.format(rootname=rootname)
            targetfile = HERRMETTARGETFILE.format(rootname=rootname)

            if not os.path.isfile(paramfile):
                raise IOError(paramfile)

            if not os.path.isfile(targetfile):
                raise IOError(paramfile)

            print(paramfile, self.paramfiles[i])
            self.paramfiles[i] = paramfile
            self.targetfiles[i] = targetfile
        self.paramfiles = np.asarray(self.paramfiles, str)
        self.targetfiles = np.asarray(self.targetfiles, str)

    def fill_extraction_files(self):

        for n, node in enumerate(self.nodes):
            rootname = os.path.join(self.rootdir, ROOTNAME.format(node=node))
            for k, p in enumerate([0.16, 0.5, 0.84]):
                f = HERRMETEXTRACTPDFMODELFILE.format(
                    rootname=rootname,
                    extract_mode=self.extract_mode,
                    extract_limit=self.extract_limit,
                    extract_llkmin=self.extract_llkmin,
                    extract_step=self.extract_step, percentile=p)
                if not os.path.isfile(f):
                    raise IOError(f)

                if k == 0:
                    self.p16files[n] = f
                elif k == 1:
                    self.medianfiles[n] = f
                elif k == 2:
                    self.p84files[n] = f

        self.medianfiles = np.asarray(self.medianfiles, str)
        self.p16files = np.asarray(self.p16files, str)
        self.p84files = np.asarray(self.p84files, str)

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def __iter__(self):
        # warning does not work for double iteration
        self._iter_index = 0
        return self

    def __next__(self):
        if self._iter_index >= len(self):
            raise StopIteration

        out = [self.nodes[self._iter_index],
               self.lons[self._iter_index],
               self.lats[self._iter_index],
               self.targetfiles[self._iter_index],
               self.paramfiles[self._iter_index],
               self.medianfiles[self._iter_index],
               self.p16files[self._iter_index],
               self.p84files[self._iter_index]]
        self._iter_index += 1
        return out

    next = __next__  # python 2

    def __len__(self):
        return len(self.nodes)

    def __str__(self):
        s = """# Nodefile for HerrMet
        # directory where inversion is run (relative to this file)
        #met rootdir = "{rootdir:s}"
        # extraction parameters needed to locate the temporary files
        #met extract_mode = "{extract_mode:s}"
        #met extract_limit = {extract_limit:d}
        #met extract_llkmin = {extract_llkmin:d}
        #met extract_step = {extract_step:d}
        #fld longitude latitude node
        #unt deg       deg      -
        #fmt %f        %f       %s""".format(
            rootdir=self.rootdir,
            extract_mode=self.extract_mode,
            extract_limit=self.extract_limit,
            extract_llkmin=self.extract_llkmin,
            extract_step=self.extract_step)
        s = "# longitude_deg latitude_deg node\n"
        for lon, lat, node in zip(self.lons, self.lats, self.nodes):
            s += "{:f} {:f} {:s}\n".format(lon, lat, node)
        return s

    def write(self, filename):
        with open(filename, 'w') as fid:
            fid.write(str(self))


class NodeFile(NodeFileString):
    def __init__(self, filename):
        """
        :param filename:
        """

        with open(filename, "r") as fid:
            NodeFileString.__init__(self,
                                    string="".join(fid.readlines()),
                                    rootdir=os.path.dirname(filename))
