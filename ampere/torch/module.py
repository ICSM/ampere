""" This module defines the base class for all `ampere.torch` objects. """

from typing import Union, List, Dict, Any
import torch
from gpytorch import Module as GPyTorchModule

class Module(GPyTorchModule):
    """
    Base class for all `ampere.torch` objects, including models, likelihoods,
    datasets and inference algorithms.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def forward(self, *args, **kwargs):
        raise NotImplementedError("The forward method must be implemented in a subclass.")

    def __call__(self, *args, **kwargs):
        return self.forward(*args, **kwargs)