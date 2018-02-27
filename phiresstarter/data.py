from __future__ import print_function

import numpy as np


class Data(object):
    """


    """

    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()


#1. Should all the photometry be stored in one object

class Photometry(Data):

    """


    """

    def __init__(self, filterName, value, uncertainty, **kwargs):
        self.filterName = filterName
        self.value = value
        self.uncertainty = uncertainty
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def synPhot(self, **kwargs):
        pass

class Spectrum(Data):

    def __init__(**kwargs):
        pass

    def __call__(self, **kwargs):
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()
    
    
