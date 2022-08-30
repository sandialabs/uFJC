"""The core module for the SWFJC single-chain model.

This module contains the core class ``SWFJC`` which, upon instantiation,
becomes a SWFJC single-chain model instance with methods for computing
single-chain quantities in either thermodynamic ensemble.
The ``SWFJCIsometric`` and ``SWFJCIsotensional`` classes are also
contained within this module.
The SWFJC is the uFJC model with a square-well link potential,
which is a special case that is efficiently treated separately
as it can be solved exactly :cite:`buche2022freely`.

Example:
    Import and create an instance of the model:

        >>> from ufjc import SWFJC
        >>> model = SWFJC()

"""

# Import internal modules
from .utility import BasicUtility

# Import external modules
import numpy as np


class SWFJCIsometric(BasicUtility):
    """The SWFJC model class for the isometric ensemble.

    Attributes:
        N (int): The number of links in the chain.
        varsigma (float): The nondimensional well width.

    """
    def __init__(self, N_b=8, varsigma=3):
        """Initializes the ``SWFJCIsometric`` class.

        Initialize and inherit all attributes and methods
        from a ``BasicUtility`` class instance.

        Args:
            N_b (int, optional, default=8):
                The number of links in the chain.
            varsigma (float, optional, default=3):
                The nondimensional well width.

        """
        super().__init__()
        self.N_b = N_b
        self.varsigma = varsigma


class SWFJCIsotensional(SWFJCIsometric):
    """The SWFJC model class for the isotensional ensemble.

    Attributes:
        N (int): The number of links in the chain.
        varsigma (float): The nondimensional well width.

    """
    def __init__(self, N_b=8, varsigma=2):
        """Initializes the ``SWFJCIsotensional`` class.

        Initialize and inherit all attributes and methods
        from a ``SWFJCIsometric`` class instance.

        Args:
            N_b (int, optional, default=8):
                The number of links in the chain.
            varsigma (float, optional, default=3):
                One plus the nondimensional well width.

        """
        super().__init__()
        self.N_b = N_b
        self.varsigma = varsigma

    def z(self, eta):
        r"""The nondimensional single-link isotensional partition function
        as a function of the nondimensional force,

        .. math::
            \mathfrak{z}(\eta) =
            \frac{1}{\eta^3}\left[
                \varsigma\eta\cosh(\varsigma\eta)
                - \sinh(\varsigma\eta)
                - \eta\cosh(\eta) + \sinh(\eta)
            \right].

        Args:
            v (array_like): The nondimensional force.

        Returns:
            numpy.ndarray: The nondimensional isotensional partition function.

        """
        eta = self.np_array(eta)
        eta_zero = eta == 0
        eta[eta_zero] = -1
        z = (
            np.sinh(eta) - eta*np.cosh(eta)
            - np.sinh(self.varsigma*eta)
            + self.varsigma*eta*np.cosh(self.varsigma*eta)
        )/eta**3
        z[eta_zero] = (self.varsigma**3 - 1)/3
        return z

    def beta_varphi(self, eta):
        r"""The nondimensional isotensional free energy as a function of
        the nondimensional force,

        .. math::
            \beta\varphi(\eta) =
            -\ln\mathfrak{z}(\eta).

        Note that this becomes the isotensional free energy
        of the FJC model as :math:`\varsigma` goes to zero,

        .. math::
            \lim_{\varsigma\to 0}\beta\varphi(\eta) =
            \ln\left[\frac{\eta}{\sinh(\eta)}\right].

        Args:
            v (array_like): The nondimensional force.

        Returns:
            numpy.ndarray: The nondimensional isotensional free energy.

        Example:
            Plot the nondimensional isotensional free energy as a function of
            the nondimensional force for varying :math:`\varsigma`:

            .. plot::

                >>> import numpy as np
                >>> import matplotlib.pyplot as plt
                >>> from ufjc.swfjc import SWFJCIsotensional
                >>> eta = np.linspace(0, 10, 1000)[1:]
                >>> _ = plt.figure()
                >>> for varsigma in [2, 3, 5, 10]:
                ...     model = SWFJCIsotensional(varsigma=varsigma)
                ...     _ = plt.plot(eta, model.beta_varphi(eta),
                ...                  label=r'$\varsigma=$'+str(varsigma))
                >>> _ = plt.plot(eta, np.log(eta/np.sinh(eta)),
                ...              'k--', label='FJC')
                >>> _ = plt.xlabel(r'$\eta$')
                >>> _ = plt.ylabel(r'$\beta\varphi$')
                >>> _ = plt.legend()
                >>> plt.show()

        """
        return (np.log(self.z(0)) - np.log(self.z(eta)))/self.varsigma

    def gamma_isotensional(self, eta):
        r"""The nondimensional end-to-end length as a function of
        the nondimensional force in the isotensional ensemble,

        .. math::
            \gamma(\eta) =
            \frac{\partial}{\partial\eta}\,\ln\mathfrak{z}(\eta) =
            \frac{
                \varsigma^2\eta\sinh(\varsigma\eta) - \eta\sinh(\eta)
            }{
                \varsigma\eta\cosh(\varsigma\eta)
                - \sinh(\varsigma\eta)
                - \eta\cosh(\eta) + \sinh(\eta)
            }

        Note that this becomes the Langevin function of the FJC model
        as :math:`\varsigma` goes to zero,

        .. math::
            \lim_{\varsigma\to 0}\gamma(\eta) =
            \coth(\eta) - \frac{1}{\eta} =
            \mathcal{L}(\eta).

        Args:
            v (array_like): The nondimensional force.

        Returns:
            numpy.ndarray: The nondimensional end-to-end length.

        Example:
            Plot the nondimensional single-chain mechanical response
            in the isotensional ensemble for varying :math:`\varsigma`:

            .. plot::

                >>> import numpy as np
                >>> import matplotlib.pyplot as plt
                >>> from ufjc.swfjc import SWFJCIsotensional
                >>> eta = np.linspace(0, 10, 1000)[1:]
                >>> _ = plt.figure()
                >>> for varsigma in [2, 3, 5, 10]:
                ...     model = SWFJCIsotensional(varsigma=varsigma)
                ...     _ = plt.plot(
                ...         model.gamma_isotensional(eta)/varsigma,
                ...         eta, label=r'$\varsigma=$'+str(varsigma))
                >>> _ = plt.plot(1/np.tanh(eta) - 1/eta, eta,
                ...              'k--', label='FJC')
                >>> _ = plt.xlabel(r'$\gamma/\varsigma$')
                >>> _ = plt.ylabel(r'$\eta$')
                >>> _ = plt.legend()
                >>> plt.show()

        """
        return (
            self.varsigma**2*eta*np.sinh(self.varsigma*eta)
            - eta*np.sinh(eta)
        )/(
            np.sinh(eta) - eta*np.cosh(eta)
            - np.sinh(self.varsigma*eta)
            + self.varsigma*eta*np.cosh(self.varsigma*eta)
        ) - 3/eta


class SWFJC(SWFJCIsotensional):
    """The SWFJC single-chain model class.

    Attributes:
        N (int): The number of links in the chain.
        varsigma (float): The nondimensional well width.

    """
    def __init__(self):
        """Initializes the ``SWFJC`` class.

        Initialize and inherit all attributes and methods
        from a ``SWFJCIsotensional`` class instance.

        Args:
            N_b (int, optional, default=8):
                The number of links in the chain.
            varsigma (float, optional, default=3):
                The nondimensional well width.

        """
        super().__init__()
