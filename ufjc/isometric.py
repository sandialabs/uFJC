"""A module for the uFJC single-chain model in the isometric ensemble.

This module consist of the class ``uFJCIsometric`` which contains
methods for computing single-chain quantities
in the isometric (constant end-to-end vector) thermodynamic ensemble.

Example:
    Import and instantiate the class:

        >>> from ufjc.isometric import uFJCIsometric
        >>> class_instance = uFJCIsometric()

"""

# Import internal modules
from .monte_carlo import MHMCMC
from .isotensional import uFJCIsotensional

# Import external modules
import numpy as np
import numpy.linalg as la


class uFJCIsometric(uFJCIsotensional):
    """The uFJC single-chain model class for the isometric ensemble.

    This class contains methods for computing single-chain quantities
    in the isometric (constant end-to-end vector) thermodynamic ensemble.
    It inherits all attributes and methods from the ``uFJCIsotensional``
    class, which inherits all attributes and methods from the
    ``BasicUtility`` class.

    """
    def __init__(self):
        """Initializes the ``uFJCIsometric`` class.

        Initialize and inherit all attributes and methods
        from a ``uFJCIsotensional`` class instance.

        """
        uFJCIsotensional.__init__(self)

    def eta_isometric(self, gamma, **kwargs):
        r"""Main function for the isometric :math:`\eta(\gamma)`.

        This is the main function utilized to compute the isometric
        nondimensional single-chain mechanical response.
        Keyword arguments specify and are passed onto the method.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to the chosen method.

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Compute the nondimensional force for an eight-link Morse-FJC at a
            nondimensional end-to-end length of 0.8 in the isometric ensemble,
            using the Legendre transformation method from the isotensional
            ensemble, and using the reduced asymptotic approach to compute
            quantities in the isotensional ensemble:

                >>> from ufjc import uFJC
                >>> model = uFJC(N_b=8, potential='morse')
                >>> model.eta_isometric([0, 0.8], \
                ...     method='legendre', approach='reduced')
                array([0.        , 4.41715473])

        Warning:
            Only the Legendre transformation method is currently unavailable:

                >>> from ufjc import uFJC
                >>> uFJC().eta_isometric(0.8, method='exact')
                array([nan])

        """
        gamma = self.np_array(gamma)
        method = kwargs.get('method', 'legendre')
        if method == 'legendre':
            return self.eta_isometric_legendre(gamma, **kwargs)
        else:
            return np.nan*gamma

    def eta_isometric_legendre(self, gamma, **kwargs):
        r"""The Legendre transformation method of approximating
        the isometric :math:`\eta(\gamma)`.

        This function uses the Legendre transformation method to obtain an
        approximate isometric nondimensional single-chain mechanical response.
        The result is to simply use the isotensional :math:`\eta(\gamma)`,
        and this approximation is asymptotically valid for :math:`N_b\gg 1`
        and appreciable loads :cite:`buche2020statistical`.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``_eta_isotensional``.

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Compute the nondimensional force at a large nondimensional
            end-to-end length using the Legendre transformation method:

                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> model.eta_isometric_legendre(1.3)
                array([28.71102552])

        """
        return self.eta_isotensional(gamma, **kwargs)

    def gamma_isometric(self, eta, **kwargs):
        r"""Main function for the isometric :math:`\gamma(\eta)`.

        This function obtains the isometric nondimensional
        single-chain mechanical response :math:`\gamma(\eta)`
        by inverting the isometric :math:`\eta(\gamma)`.

        Args:
            eta (array_like): the nondimensional force(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``_eta_isometric``.

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Check that :math:`\gamma[\eta(\gamma)] = \gamma\,`:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> def check_eta(gamma):
                ...     eta_fun = lambda gamma: model.eta_isometric(gamma)
                ...     gamma_fun = lambda eta: model.gamma_isometric(eta)
                ...     return np.isclose(gamma_fun(eta_fun(gamma))[0], gamma)
                >>> check_eta(np.random.rand())
                True

        """
        def eta_fun(gamma):
            return self.eta_isometric(gamma, **kwargs)
        return self.inv_fun_1D(eta, eta_fun)

    def vartheta_isometric(self, gamma, **kwargs):
        r"""Main function for the isometric :math:`\vartheta(\gamma)`.

        This is the main function utilized to compute the nondimensional
        Helmholtz free energy per link, an isometric quantity.
        Keyword arguments specify and are passed onto the method.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to the chosen method.

        Returns:
            numpy.ndarray: The nondimensional Helmholtz free energy per link.

        Example:
            Compute the nondimensional Helmholtz free energy per link
            for an eight-link Morse-FJC at a
            nondimensional end-to-end length of 0.8 in the isometric ensemble,
            using the Legendre transformation method from the isotensional
            ensemble, and using the reduced asymptotic approach to compute
            quantities in the isotensional ensemble:

                >>> from ufjc import uFJC
                >>> model = uFJC(N_b=8, potential='morse')
                >>> model.vartheta_isometric(0.8, \
                ...     method='legendre', approach='reduced')
                array([1.23847534])

        Warning:
            The exact method is currently unavailable:

                >>> from ufjc import uFJC
                >>> uFJC().vartheta_isometric(0.8, method='exact')
                nan

        """
        method = kwargs.get('method', 'legendre')
        if method == 'exact':
            return np.nan*gamma
        elif method == 'legendre':
            return self.vartheta_isometric_legendre(gamma, **kwargs)

    def vartheta_isometric_legendre(self, gamma, **kwargs):
        r"""The Legendre transformation method of approximating
        the isometric :math:`\vartheta(\gamma)`.

        This function uses the Legendre transformation method to obtain an
        approximate isometric Helmholtz free energy per link.
        The result is independent of the number of links :math:`N_b`, and
        this approximation is asymptotically valid for :math:`N_b\gg 1`
        and appreciable loads :cite:`buche2021chain`.
        For example, using the reduced asymptotic approach, this is

        .. math::
            \vartheta(\gamma) \sim
            \ln\left\{\frac{
            \eta\exp[\eta\mathcal{L}(\eta)]}{\sinh(\eta)}\right\}
            + \beta u[\lambda(\eta)]
            ,

        valid when :math:`\varepsilon\gg 1` and :math:`N_b\gg 1` are
        simultaneously true.
        Note that :math:`\eta=\eta(\gamma)` is implied, and obtained
        through inverting the isotensional :math:`\gamma(\eta)`.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to the chosen method.

        Returns:
            numpy.ndarray: The nondimensional Helmholtz free energy per link.

        Example:
            Approximate the nondimensional Helmholtz free energy per link
            using the Legendre method and both asymptotic approaches:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='log-squared', varepsilon=23)
                >>> model.vartheta_isometric_legendre(1.1)
                array([1.90431381])
                >>> model.vartheta_isometric_legendre(1.1, approach='reduced')
                array([2.09238198])

        Warning:
            Only the asymptotic approaches are currently unavailable:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='log-squared', varepsilon=23)
                >>> model.vartheta_isometric_legendre(1.1, approach='exact')
                nan

        """
        # Invert gamma=gamma(eta) for the corresponding eta
        eta = self.eta_isotensional(gamma, **kwargs)

        # Avoid overflow, important for integrating P_eq
        eta[eta > self.maximum_exponent] = self.maximum_exponent

        # Find the corresponding bond stretch under direct eta
        lambda_ = 1 + self.delta_lambda(eta)

        # Need to finish this portion
        approach = kwargs.get('approach', 'asymptotic')
        if approach == 'asymptotic':
            Ln = self.langevin(eta)
            coth = self.coth(eta)
            return eta*Ln + self.log_over_sinh(eta) + \
                self.varepsilon*(self.phi(lambda_) - self.phi(1)) + \
                eta**2/self.kappa*(
                    (1 - Ln*coth)/(self.c + eta/self.kappa*coth)
                ) - np.log(1 + eta*coth/self.kappa)
        elif approach == 'reduced':
            return eta*self.langevin(eta) + self.log_over_sinh(eta) + \
                self.varepsilon*(self.phi(lambda_) - self.phi(1))
        else:
            return np.nan*gamma

    def beta_U_config(self, config):
        r"""The nondimensional potential energy of a configuration.

        This function provides the nondimensional potential energy
        :math:`\beta U` given the configuration of the chain, i.e. the
        vector position of each atom/hinge relative to the first one.

        Args:
            config (numpy.ndarray): The configuration of the chain,
                a :math:`(N_b+1)`-by-3 numpy array.

        Returns:
            float: The nondimensional potential energy :math:`\beta U`.

        Example:
            Compute the potential energy of the uniformly-stretched
            default initial configuration:

                >>> from ufjc import uFJC
                >>> model = uFJC(N_b=8, potential='lennard-jones')
                >>> model.beta_U_config(1.1*model.init_config)
                -570.4631978476273

        """
        beta_U = 0
        for j in range(1, len(config)):
            lambda_ = la.norm(config[j, :] - config[j - 1, :])
            beta_U += self.beta_u(lambda_)
        return beta_U
