class Relation(object):
    # I need a pickable object to pass a function that was defined from a string...
    def __init__(self, name, string):
        self.name = name
        self.string = string
        assert self.string.strip().startswith('def {}('.format(name))
        self.fun = None

    def __getstate__(self):
        return self.name, self.string

    def __setstate__(self, state):
        name, string = state
        Relation.__init__(self, name=name, string=string)

    def __call__(self, *args, **kwargs):
        try:
            return self.fun(*args, **kwargs)
        except TypeError as e:
            if str(e) == "'NoneType' object is not callable":

                import imp
                modulename = "{}".format(self.name)
                pyfilename = "./_relation_{}.py".format(modulename)
                funcname = self.name

                with open(pyfilename, 'w') as fid:
                    fid.write(self.string + "\n")

                self.fun = getattr(imp.load_source(modulename, pyfilename), funcname)
                return self.fun(*args, **kwargs)

            else:
                raise e