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

    
    """
    These are the arithmetic operations. These should accept two kinds of "other" arguments 
    
    1) Other Model (or subclass) instances
    2) python/numpy numeric types (i.e. all types of floats and ints)

    and everything else should return NotImplemented. In these two cases, the methods should
    return a new CompositeModel instance that will provide the interface to the pair of models 
    and combine their ModelFluxes with the appropriate operation. The CompositeModel instance
    is repsonsible for making sure that the wavelengths of the two Models are the same. More 
    complexity may be required when we have spatial information to worry about as well.

    def __mul__(self, other):
        pass
    """
    def __add__(self, other): #This method will be called if x + y is invoked and x is a Model instance  (self + other)
        if isinstance(other, Model): #When we try to combine two models we need to identify it properly
            newModel = CompositeModel(self, other, '+')
            pass
        else:
            return NotImplemented

    def __radd__(self, other):  #This method will be called if x + y is invoked and x is not a Model instance and x.__add__ can't handle being given a Model instance (other + self)
        return NotImplemented

    def __iadd__(self, other): #This method will be called if x+=y is called and x is a Model instance (self += other)
        return NotImplemented

    def __sub__(self, other): #This method is like __add__ but for x - y
        return NotImplemented

    def __rsub__(self, other): #This method is like __radd__ but for x - y
        return NotImplemented

    def __isub__(self, other): #This method is like __iadd__ but for x -= y
        return NotImplemented
    
    def __mul__(self, other): #This method is like __add__ but for x * y
        return NotImplemented

    def __rmul__(self, other): #This method is like __radd__ but for x * y
        return NotImplemented

    def __imul__(self, other): #This method is like __iadd__ but for x *= y
        return NotImplemented

    def __truediv__(self, other): #This method is like __add__ but for x / y
        return NotImplemented

    def __rtruediv__(self, other): #This method is like __radd__ but for x / y
        return NotImplemented

    def __itruediv__(self, other): #This method is like __iadd__ but for x /= y
        return NotImplemented

    def lnprior(self, **kwargs):
        raise NotImplementedError()

class CompositeModel(Model):
    ''' 
    Class which represents models which are composed of arithmetic combinations of other models

    The basic idea is that this class simply handles the bookkeeping of combining the results of two
    Model instances with a given operator. Both must be defined on the same wavelengths and model
    the same kind of observable (i.e. for now, spectral flux density as a function of wavelength, 
    not images). 
    '''

    def __init__(self, model1, model2, operator, **kwargs):
        self.model1 = model1
        self.model2 = model2
        self.operator = operator
        #Now look up how many arguments each one takes and define the same things that all models need to define
        self.npars = model1.npars + model2.npars #how many arguments do we have to accept in total
        self.part1 = model1.npars #How many of those arguments go to model 1


        self.wavelength = self.model1.wavelength
        raise NotImplementedError()

    def __call__(self, *args, **kwargs): #This method should invoke the __call__ routines of the two underlying Model instances and combine their ModelFluxes with the specified operator
        self.model1(args[:self.part1])
        self.model2(args[self.part1:])
        self.modelFlux = eval("self.model1.modelFlux"+self.operator+"self.model1.modelFlux") #this needs to be generalised to models that have more complex predictions than just a spectrum
        #raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def lnprior(self, theta, **kwargs): #This method should simply call the priors of the underlying Models for the appropriate chunks of theta and return their sum
        return self.model1.lnprior(theta[:self.part1]) + self.model2.lnprior(theta[self.part1:])
    
    def prior_transform(self, u, **kwargs): #This should split u into two components of the size required by each Model and feed them in, then reconstruct the required format of the resulting parameter values
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
