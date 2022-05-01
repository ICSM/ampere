from __future__ import print_function

import numpy as np
from inspect import signature
from .basesearch import BaseSearch


class BaseNestedSampler(Basesearch):

    def __init__(self, model = None, data= None, **kwargs):
        self.model = model
        self.dataSet = data
