import os
import numpy as np
from srfpython.HerrMet.files import ROOTNAME, HERRMETPARAMFILE, HERRMETTARGETFILE, HERRMETEXTRACTPDFMODELFILE


"""
a node file is an ascii file with 3 columns and one line header 

# longitude_deg latitude_deg node
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
        assert os.path.isdir(rootdir)
        self.rootdir = rootdir
        lines = string.split('\n')

        lons, lats, nodes = [], [], []
        paramfiles = []
        targetfiles = []

        for line in lines:
            line = line.strip()
            if line.startswith("#"):
                continue

            if not len(line):
                continue

            lon, lat, node = line.split()
            lons.append(lon)
            lats.append(lat)
            nodes.append(node)

            rootname = os.path.join(self.rootdir, ROOTNAME.format(node=node))
            paramfile = HERRMETPARAMFILE.format(rootname=rootname)
            targetfile = HERRMETTARGETFILE.format(rootname=rootname)

            if not os.path.isfile(paramfile):
                raise IOError(paramfile)

            if not os.path.isfile(targetfile):
                raise IOError(paramfile)

            paramfiles.append(paramfile)
            targetfiles.append(targetfile)

        self.lons = np.asarray(lons, float)
        self.lats = np.asarray(lats, float)
        self.nodes = np.asarray(nodes, str)
        self.paramfiles = np.asarray(paramfiles, str)
        self.targetfiles = np.asarray(targetfiles, str)
        self.medianfiles = np.zeros(len(lons), object)  # not str otherwise => type is set to |S1
        self.p16files = np.zeros(len(lons), object)
        self.p84files = np.zeros(len(lons), object)

        if not len(self.nodes) == len(np.unique(self.nodes)):
            raise ValueError('node names are not unique')

    def fill_extraction_files(
            self, extract_mode="best", extract_limit=1000,
            extract_llkmin=0, extract_step=1):
        """
        :param extract_mode:
        :param extract_limit:
        :param extract_llkmin:
        :param extract_step:
        :return:
        """

        for n, node in enumerate(self.nodes):
            rootname = os.path.join(self.rootdir, ROOTNAME.format(node=node))
            for k, p in enumerate([0.16, 0.5, 0.84]):
                f = HERRMETEXTRACTPDFMODELFILE.format(
                    rootname=rootname, extract_mode=extract_mode,
                    extract_limit=extract_limit, extract_llkmin=extract_llkmin,
                    extract_step=extract_step, percentile=p)
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

    def __iter__(self):
        self._current = 0
        return self

    def __next__(self):
        if self._current >= len(self):
            raise StopIteration

        out = [self.nodes[self._current],
               self.lons[self._current],
               self.lats[self._current],
               self.targetfiles[self._current],
               self.paramfiles[self._current],
               self.medianfiles[self._current],
               self.p16files[self._current],
               self.p84files[self._current]]
        self._current += 1
        return out

    next = __next__  # python 2

    def __len__(self):
        return len(self.nodes)

    def __str__(self):
        s = "# longitude_deg latitude_deg node\n"
        for lon, lat, node in zip(self.lons, self.lats, self.nodes):
            s += "{:f} {:f} {:s}\n".format(lon, lat, node)
        return s

    def write(self, filename):
        with open(filename, 'w') as fid:
            fid.write(str(self))


class NodeFile(NodeFileString):
    def __init__(self, filename, rootdir=None):
        """
        :param filename:
        """
        if rootdir is None:
            rootdir = os.path.dirname(filename)

        with open(filename, "r") as fid:
            NodeFileString.__init__(self, string="".join(fid.readlines()), rootdir=rootdir)
