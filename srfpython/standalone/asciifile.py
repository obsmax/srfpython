from __future__ import print_function

import numpy as np
from numpy import array, zeros, ones #needed for eval

"""
an object to manage ascii files with data, metadata with different types
"""


demo="""
#demonstration file
#comments order will not be preserved 
#indicate metadata with key word #met

#met a="array([1])"  
#met b=array([1])    
#met c=124.2
#met d=[1,2,3]

#indicate the field names, units and formating using #fld, #unt and #fmt respectively
#fld FIELD1 FIELD2 FIELD3
#unt km     -      -
#fmt %.0f   %1s    %10s
     1.     A      lklkj      
     2.     B      123
     
     3.     C      kqklj
"""


class AsciiFile_fromstring(object):

    def __init__(self, string):
        self.metadata = {} #dictionnary with metadata evaluated from ascii file
        self.data = None #structured array, see https://docs.scipy.org/doc/numpy-1.13.0/user/basics.rec.html
        self.fields = None
        self.formats = None
        self.fmtline = ""
        self.units = None
        self.comments = []
        dtype = None
        for l in string.split('\n'):
            l = l.strip()
            if l == "": continue
            elif l.startswith("#"):
                if l.startswith("#met"):
                    l = l.lstrip('#met').split('#')[0]
                    #key, value = [_.strip() for _ in l.split('=')]
                    key = l.split('=')[0].strip()
                    value = "=".join(l.split('=')[1:])
                    self.metadata[key] = eval(value)
                elif l.startswith('#fld'):
                    self.fields = np.asarray(l.lstrip("#fld").split(), str)
                    self.M = len(self.fields)
                elif l.startswith('#unt'):
                    self.units = np.asarray(l.lstrip("#unt").split(), str)
                elif l.startswith('#fmt'):
                    self.formats = np.asarray(l.lstrip("#fmt").split(), str)
                    self.fmtline = l.replace('#fmt', "    ")
                    self.types = np.empty(len(self.formats), object)
                    for n, fmt in enumerate(self.formats):
                        if fmt.endswith('s'):
                            self.types[n] = "S100"
                        elif fmt.endswith('f'):
                            self.types[n] = float
                        elif fmt.endswith('d'):
                            self.types[n] = int
                        else:
                            raise Exception('unknown column format %s' % fmt)
                else:
                    self.comments.append(l)

            else:
                assert self.fields is not None
                assert self.formats is not None
                dtype = [_ for _ in zip(self.fields, self.types)]
                values = tuple(l.split())
                assert len(values) == self.M
                if self.data is None:
                    self.data = [values]
                else:
                    self.data.append(values)
        self.data = np.array(self.data, dtype=dtype)

    def __len__(self):
        return len(self.data)

    def __str__(self):
        out = ""
        if len(self.comments):
            out += "\n".join(self.comments) + "\n"
        for k, v in self.metadata.items():
            if isinstance(v, str):
                out += "#met %s = '%s'\n" % (k, str(v))
            else:
                out += "#met %s = %s\n" % (k, str(v))
        out += "#fld " + " ".join(self.fields) + "\n"
        out += "#unt " + " ".join(self.units) + "\n"
        out += "#fmt " + self.fmtline.lstrip() + "\n"
        for i, line in enumerate(self.data):
            out += self.fmtline % tuple(line) + "\n"
        return out

    def write(self, filename=None):
        if filename is None:
            print (self.__str__())
        else:
            with open(filename, 'w') as fid:
                fid.write(self.__str__())

    def __getitem__(self, item):
        if isinstance(item, str):
            "item is a field name, return the column"
            return self.data[item]
        elif isinstance(item, int):
            "item is a row number, return the line"
            return self.data[item]
        raise IndexError('indexing with object of type %s is not implemented' % str(type(item)))


class AsciiFile(AsciiFile_fromstring):
    def __init__(self, filename):
        with open(filename, 'r') as fid:
            string = "".join(fid.readlines())
        AsciiFile_fromstring.__init__(self, string)


if __name__ == "__main__":
    a = AsciiFile_fromstring(demo)
    print(a.metadata)
    print ("---")
    print (a.data)
    print ("---")
    print (a)

