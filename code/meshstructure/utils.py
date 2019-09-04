class lazyattr(object):
    """A read-only attribute that is created on-demand. Like @property
    but only evaluates the body of the function once."""
    def __init__(self, fget, doc=None):
        self.fget = fget
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__
        self.__module__ = fget.__module__

    def __get__(self, obj, cls):
        if obj is None:
            return self
        return obj.__dict__.setdefault(self.__name__, self.fget(obj))
