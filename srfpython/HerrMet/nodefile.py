import numpy as np
import os

"""
a node file is an ascii file with 3 columns and one line header 

# longitude_deg latitude_deg node
12.132131 45.355425434 node0001
22.132131 35.355425434 node0002
32.132131 25.355425434 node0003
...
"""


class NodeFileString(object):
    def __init__(self, string):
        """

        :param string:
        """

        lines = string.split('\n')

        lons, lats, nodes = [], [], []
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

        self.lons = np.asarray(lons, float)
        self.lats = np.asarray(lats, float)
        self.nodes = np.asarray(nodes, str)

        assert len(self.nodes) == len(np.unique(self.nodes))

    def __str__(self):
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
            NodeFileString.__init__(self, "".join(fid.readlines()))

