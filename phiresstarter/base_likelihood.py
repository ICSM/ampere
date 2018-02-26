from __future__ import print_function

import numpy as np



class BaseLikelihood(object):
    ''' docstring goes here  '''

    def __init__(self, theta, **kwargs):
        self.theta = theta
        pass

    def __repr__(self, **kwargs):
        return 'Shaddapayouface'

    def __str__(self, **kwargs):
        return 'What do you want?'

    def lnprior(self, theta, **kwargs):
        a = theta
        return 0

    def lnlike(self, theta, **kwargs):
        return 0

    def lnprob(self, theta, **kwargs):
        lnprob = self.lnprior(theta, **kwargs)
        return lnprob + self.lnlike(theta, **kwargs)



if __name__ == "__main__":
    a = BaseLikelihood(1)
    b = BaseLikelihood(2)
    print(a.theta)
    #help(a)
    dir(a)
    theta = 1
    print(a.lnprob(theta))
