class Relation(object):
    # I need a pickable object to pass a function that was defined from a string...
    def __init__(self, name, string):
        """
        :param name: name of the relation (i.e. function name)
        :param string: expression of the function
            version 1 : string = "lambda x: 3 * x + 2"
            version 2 : string = "def name(x): return 3 * x + 2"
        """
        self.name = name
        self.string = string

        if self.string.strip().startswith('lambda'):
            self.fun = eval(self.string)

        elif self.string.strip().startswith('def {}('.format(name)):
            assert "return" in self.string
            self.fun = eval(compile(self.string, '<string>', 'exec'))

        else:
            raise ValueError("could not compile {} ".format(string))

    def __getstate__(self):
        return self.name, self.string

    def __setstate__(self, state):
        name, string = state
        Relation.__init__(self, name=name, string=string)

    def __call__(self, *args, **kwargs):
        return self.fun(*args, **kwargs)
