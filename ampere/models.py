import numpy as np


class Model(object):
    '''
    Base class on which to build different kinds of models.

    '''

    def __init__(self, **kwargs):
        raise NotImplementedError()

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()




class AnalyticalModel(Model):
    '''
    Class on which to build analytical models.

    '''

    def __init__(self, **kwargs):
        raise NotImplementedError()

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()


class RTModel(Model):
    '''
    Class on which to build radiative-transfer models.

    '''

    def __init__(self, **kwargs):
        raise NotImplementedError()

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def inputFile(self, **kwargs):
        raise NotImplementedError()
