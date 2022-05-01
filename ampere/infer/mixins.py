""" Mixins for inference

This module contains a set of mixins used by Ampere's inference classes to 
enable smoother composition of samplers. The idea is to maximise the sharing of 
code without ending up with too many layers of single inheretance and avoiding 
problems with multiple inheretance. 

Classes defined here *MUST* be defined as mixins only. That means their only 
parent class should be object (or another mixin whose methods they overload) and
 they must *NOT* define an __init__() method, or any methods defined in 
BaseSearch or any of its children. All class names should have "Mixin" appended
to them so it is obvious what they are being used for when defining other 
classes
"""

from __future__ import print_function

import numpy as np




