"""

"""

# Import external modules
import numpy as np


class SWFJCIsometric(object):
    """

    """
    def __init__(self):
        """

        """
        pass


class SWFJCIsotensional(SWFJCIsometric):
    """

    """
    def __init__(self, N_b=8, varsigma=3):
        """

        """
        super().__init__()
        self.N_b = N_b
        self.varsigma = varsigma

    def z(self, eta):
        """

        """
        return (
            np.sinh(eta) - eta*np.cosh(eta)
            - np.sinh((1 + self.varsigma)*eta)
            + (1 + self.varsigma)*eta*np.cosh((1 + self.varsigma)*eta)
        )/eta**3

    def varphi(self, eta):
        """

        """
        return np.log(self.z(eta)) - np.log(self.z(0))

    def gamma_isotensional(self, eta):
        r"""The nondimensional end-to-end length as a function of
        the nondimensional force in the isotensional ensemble,

        .. math::
            \gamma(\eta) =
            -\frac{\partial}{\partial\eta}\,\beta\varphi(\eta).

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
                >>> from swfjc import SWFJCIsotensional
                >>> eta = np.linspace(0, 10, 1000)[1:]
                >>> _ = plt.figure()
                >>> for varsigma in [0.001, 1, 3, 10, 30]:
                ...     model = SWFJCIsotensional(varsigma=varsigma)
                ...     _ = plt.plot(model.gamma_isotensional(eta), eta,
                ...                  label=r'$\varsigma=$'+str(varsigma))
                >>> _ = plt.plot(1/np.tanh(eta) - 1/eta, eta,
                ...              'k--', label='FJC')
                >>> _ = plt.xlabel(r'$\gamma$')
                >>> _ = plt.ylabel(r'$\eta$')
                >>> _ = plt.legend()
                >>> plt.show()

        """
        return (
            3*np.sinh(eta)
            - 3*np.sinh((1 + self.varsigma)*eta)
            - 3*eta*np.cosh(eta)
            + 3*(1 + self.varsigma)*eta*np.cosh((1 + self.varsigma)*eta)
            + eta**2*np.sinh(eta)
            - ((1 + self.varsigma)*eta)**2*np.sinh((1 + self.varsigma)*eta)
        )/(
            (1 + self.varsigma)*eta*(
                np.sinh((1 + self.varsigma)*eta)
                - (1 + self.varsigma)*eta*np.cosh((1 + self.varsigma)*eta)
                - np.sinh(eta)
                + eta*np.cosh(eta)
            )
        )


class SWFJC(SWFJCIsotensional):
    """

    """
    def __init__(self):
        """

        """
        super().__init__()

    def gamma(self, eta, ensemble='isotensional'):
        """
        
        """
        if ensemble == 'isotensional':
            return self.gamma_isotensional(eta)
        elif ensemble == 'isometric':
            return self.gamma_isometric(eta)

    def eta(self, gamma, ensemble='isometric'):
        """
        
        """
        if ensemble == 'isotensional':
            return self.eta_isotensional(gamma)
        elif ensemble == 'isometric':
            return self.eta_isometric(gamma)
