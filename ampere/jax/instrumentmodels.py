""" This module defines the base class for instrument models in `ampere.jax`,
along with a few example instrument models. """ 

from typing import Union, List, Dict, Any
import jax
import jax.numpy as jnp
from .module import Module
from .transformations import Transformation

class InstrumentModel(Module):
    """
    Base class for all instrument models in `ampere.jax`. This class is intended
    to be used as a base class for all instrument models, and provides a
    `__call__` method that should be implemented in a subclass.
    """

    transformations: List[Transformation]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError("The __call__ method must be implemented in a subclass.")
    
    def add_transformation(self, transformation: Transformation):
        self.transformations.append(transformation)