import .models

class Star(AnalysticalModel):
    pass



class Tlusty(Star):
    '''
    Base class of a TLUSTY stellar model. 
    '''
    def __init__(self, **kwargs):
        pass


    
class DustyTlusty(Tlusty):
    '''
    Class for modelling Star + extinction

    Assumes that you know the spectral type (and hence temperature) and want to model the extinction curve to constrain the properties of the dust.
    '''

    def __init__(self, starFile,flatprior=True,
                 **kwargs):
        """ First read the model spectrum and do the unit conversion """


        """  """
        pass

    def __call__(self, theta, **kwargs):
        """ Compute dust model """


        """ Extinguish basic stellar model with dust model """

        
        raise NotImplementedError()
    
    def __str__(self, **kwargs):
        raise NotImplementedError()

    def __repr__(self, **kwargs):
        raise NotImplementedError()
