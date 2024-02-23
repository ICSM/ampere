""" This module defines the base class for all `ampere.jax` models. """

from typing import Union, List, Dict, Any
import jax
import jax.numpy as jnp
from jax import vmap
from jaxtyping import Array, PyTree, Float
from functools import partial
import equinox as eqx
from equinox import Module as EqxModule
from equinox import AbstractVar
from equinox import Dataset as EqxDataset
from equinox import Parameter as EqxParameter


class Module(EqxModule):
    """
    Base class for all `ampere.jax` objects, including models, likelihoods,
    datasets and inference algorithms.
    """

    # We need to define some abstract variables to make sure that the `Module` class 
    # is properly implemented.
    parameters: AbstractVar[Union[List[EqxParameter], Dict[str, EqxParameter]]]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def forward(self, *args, **kwargs):
        raise NotImplementedError("The forward method must be implemented in a subclass.")

    def __call__(self, *args, **kwargs):
        return self.forward(*args, **kwargs)
    

