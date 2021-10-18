try:
    from ._cmakeinfo import ncplugin_name as _name
except ImportError:
    raise ImportError("Could not import . _cmakeinfo (you should first install with CMake before running)")
from .cbindings import hooks as _hooks
import numpy as np
import random
random.seed(123456)#reproducible random streams

class PhysicsModel:

    #"""Python version of the PhysicsModel class defined in the plugin."""
    #def __init__(self,A, s, rg, m, p, sigma0):
    #    """Initialise model. Refer to NCPhysicsModel.hh for the meaning of the parameters."""
    #    self.__modelparams=(A, s, rg, m, p, sigma0)

    def __init__(self,filename):
        """Initialise model. Refer to NCPhysicsModel.hh for the meaning of the parameters."""
        self.__modelparams=filename
    def calcCrossSection(self,ekin):
        """Calculates cross sections. The ekin parameter can be a float or a numpy
        array, and the return value will be similar."""
        scalar = not hasattr(ekin,'__len__')
        ekin = np.atleast_1d(ekin)
        xs = np.empty(ekin.size)
        _hooks['getmanyxsvalues'](*self.__modelparams,ekin.size,ekin,xs)
        return xs[0] if scalar else xs

    def sampleScatMu(self,ekin,nvalues = 1):
        """Samples scattering mu=cos(theta_scat). The ekin parameter must be a float. If
        nvalues>1, returns a numpy array of results."""
        mu = np.empty(nvalues)
        _hooks['samplemanyscatmu'](*self.__modelparams,ekin,nvalues,mu)
        return mu if nvalues>1 else mu[0]

def registerPlugin():
    """Activate the plugin by registering with NCrystal."""
    import NCrystal#Make sure NCrystal is loaded first
    _hooks['ncplugin_register']()

def pluginName():
    return _name

