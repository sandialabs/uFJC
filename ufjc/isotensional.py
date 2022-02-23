"""A module for the uFJC single-chain model in the isotensional ensemble.

This module consist of the class ``uFJCIsotensional`` which contains
methods for computing single-chain quantities
in the isotensional (constant end-to-end force) thermodynamic ensemble.

Example:
    Import and instantiate the class:

        >>> from ufjc.isotensional import uFJCIsotensional
        >>> class_instance = uFJCIsotensional()

"""

# Import internal modules
from .monte_carlo import MHMCMC
from .utility import BasicUtility

# Import external modules
import numpy as np
import numpy.linalg as la
from scipy.special import erf
from scipy.integrate import quad_vec


class uFJCIsotensional(BasicUtility):
    """The uFJC single-chain model class for the isotensional ensemble.

    This class contains methods for computing single-chain quantities
    in the isotensional (constant end-to-end force) thermodynamic ensemble.
    It inherits all attributes and methods from the ``BasicUtility`` class.

    """
    def __init__(self):
        """Initializes the ``uFJCIsotensional`` class.

        Initialize and inherit all attributes and methods
        from a ``BasicUtility`` class instance.

        """
        BasicUtility.__init__(self)

    def gamma_isotensional(self, eta, **kwargs):
        r"""Main function for the isotensional :math:`\gamma(\eta)`.

        This is the main function utilized to compute the isotensional
        nondimensional single-chain mechanical response.
        Keyword arguments specify and are passed onto the approach.

        Args:
            eta (array_like): The nondimensional force(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``gamma_isotensional_MHMCMC`` when the
                Monte Carlo approach is chosen.

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Compute the nondimensional end-to-end length of the EFJC
            in the isotensional ensemble under a nondimensional force of 23
            using numerical quadrature:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='harmonic')
                >>> model.gamma_isotensional(23, approach='quadrature')
                array([1.22689438])

        Note:
            An improperly-specified approach will result in nan:

                >>> from ufjc import uFJC
                >>> uFJC().gamma_isotensional(0.8, approach='blah blah')
                array([nan])

        """
        eta = self.np_array(eta)
        approach = kwargs.get('approach', 'asymptotic')
        if approach == 'asymptotic':
            return self.gamma_isotensional_asymptotic(eta)
        elif approach == 'reduced':
            return self.gamma_isotensional_asymptotic_reduced(eta)
        elif approach == 'exact':
            return self.gamma_isotensional_exact(eta)
        elif approach == 'quadrature':
            return self.gamma_isotensional_quadrature(eta)
        elif approach == 'monte carlo':
            return self.gamma_isotensional_MHMCMC(eta, **kwargs)
        else:
            return np.nan*eta

    def gamma_isotensional_asymptotic(self, eta):
        r"""The full asymptotic approach for the
        isotensional :math:`\gamma(\eta)`.

        This function provides the asymptotic approach for
        obtaining :math:`\gamma(\eta)` in the isotensional ensemble
        :cite:`buche2022on`,
        a composite approximation for all :math:`\eta`
        that is valid for :math:`\varepsilon\gg 1`,

        .. math::
            \gamma(\eta) \sim \mathcal{L}(\eta) +
                \frac{\eta}{\kappa}\left[
                \frac{1 - \mathcal{L}(\eta)\coth(\eta)}
                {c + (\eta/\kappa)\coth(\eta)}\right]
                + \Delta\lambda(\eta)
            .

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Compute the asymptotically-correct nondimensional end-to-end length
            of the Mie-FJC in the isotensional ensemble
            under a nondimensional force of 23:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='mie')
                >>> model.gamma(23, approach='asymptotic')
                array([0.96204056])

        """
        Ln = self.langevin(eta)
        coth = self.coth(eta)
        return Ln + self.delta_lambda(eta) + \
            eta/self.kappa*((1 - Ln*coth)/(self.c + eta/self.kappa*coth))

    def gamma_isotensional_asymptotic_reduced(self, eta):
        r"""The reduced asymptotic approach for the
        isotensional :math:`\gamma(\eta)`.

        This function provides the reduced asymptotic approach for
        obtaining :math:`\gamma(\eta)` in the isotensional ensemble,
        a composite approximation for all :math:`\eta`
        that is valid for :math:`\varepsilon\gg 1`,

        .. math::
            \gamma(\eta) \sim \mathcal{L}(\eta) + \Delta\lambda(\eta).

        The simplification involves neglecting terms that contribute negligibly
        apart when :math:`\eta=\mathrm{ord}(\varepsilon)`, a region that tends
        to vanish for :math:`\varepsilon\gg 1`.
        This simplification is then asymptotically valid in the same limit that
        made the original asymptotic approximation possible.
        Further, this asymptotic approximation is equivalent to asymptotically
        matching the low-force entropic mechanical response
        :math:`\mathcal{L}(\eta)` to the high-force enthalpic
        mechanical response :math:`\lambda(\eta)` :cite:`buche2021chain`.

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Verify that either asymptotic approximation obtains close-to-exact
            end-to-end lengths over many nondimensional forces
            when the nondimensional energy scale is high:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC(potential='morse', varepsilon=888)
                >>> eta_trial = np.logspace(-2, 1, 88)
                >>> gamma_quad = model.gamma(eta_trial, approach='quadrature')
                >>> gamma_asymp = model.gamma(eta_trial, approach='asymptotic')
                >>> (np.abs(gamma_quad-gamma_asymp)/gamma_quad < 1e-2).all()
                True
                >>> gamma_simpl = model.gamma(eta_trial, approach='asymptotic')
                >>> (np.abs(gamma_quad-gamma_simpl)/gamma_quad < 1e-2).all()
                True

        """
        return self.langevin(eta) + self.delta_lambda(eta)

    def gamma_isotensional_exact(self, eta):
        r"""Exact approach for the
        isotensional :math:`\gamma(\eta)`.

        This function provides the exact, analytical approach for obtaining
        :math:`\gamma(\eta)` in the isotensional ensemble.

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Warning:
            This is currently only available for harmonic link potentials,
            and will therefore return nan for non-harmonic potentials:

                >>> from ufjc import uFJC
                >>> uFJC().gamma([88e-2, 8], approach='exact')
                array([0.29518908, 0.97632595])
                >>> uFJC(potential='morse').gamma([88e-2, 8], approach='exact')
                array([nan, nan])

        """
        if self.potential == 'harmonic':
            return self.efjc_gamma_isotensional_exact(eta)
        else:
            return np.nan*eta

    def efjc_gamma_isotensional_exact(self, eta):
        r"""The exact approach for the EFJC
        isotensional :math:`\gamma(\eta)`.

        This function provides the exact, analytical approach for obtaining
        :math:`\gamma(\eta)` in the isotensional ensemble
        for harmonic link potentials, i.e. for the EFJC model,

        .. math::
            \gamma(\eta) =
            \mathcal{L}(\eta) +
            \frac{\eta}{\kappa}\left[1 +
            \frac{1 - \mathcal{L}(\eta)\coth(\eta)}
            {1 + (\eta/\kappa)\coth(\eta)}\right]
            + \frac{\partial}{\partial\eta}\,\ln\left[1+g(\eta)\right]
            ,

        where the function :math:`g(\eta)` is given by :cite:`buche2022on`

        .. math::
            g(\eta) =
            \frac{
            e^{\eta}\left(\frac{\eta}{\kappa} + 1\right)
                \,\mathrm{erf}\left(\frac{\eta+\kappa}{\sqrt{2\kappa}}\right)
            - e^{-\eta}\left(\frac{\eta}{\kappa} - 1\right)
                \,\mathrm{erf}\left(\frac{\eta-\kappa}{\sqrt{2\kappa}}\right)}
            {4\sinh(\eta)\left[1 + (\eta/\kappa)\coth(\eta)\right]}
            - \frac{1}{2}
            .

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Compute the exact nondimensional end-to-end length of the EFJC
            in the isotensional ensemble under a nondimensional force of 23:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='harmonic')
                >>> model.gamma(23, approach='exact')
                array([1.22689438])

        """
        # Avoid overflow
        eta[eta > self.maximum_exponent] = self.maximum_exponent
        eta[eta == 0] = 1e-88
        eta[eta < 1e-88] = 1e-88

        # Pre-compute some common terms
        ep = erf((eta + self.kappa)/np.sqrt(2*self.kappa))
        en = erf((eta - self.kappa)/np.sqrt(2*self.kappa))
        tp = np.exp(eta)*(eta/self.kappa + 1)*ep
        tn = np.exp(-eta)*(eta/self.kappa - 1)*en
        tc = tp - tn
        coth = self.coth(eta)
        csch = coth/np.cosh(eta)

        # Compute the numerator and denominator of the complicated term
        denominator = 1 + eta/self.kappa*coth + tc/np.sinh(eta)/2
        numerator = (coth - eta*csch**2)/self.kappa
        numerator += -coth*csch*tc/2
        numerator += csch*np.exp(eta)*(eta/self.kappa + 1 + 1/self.kappa)*ep/2
        numerator += csch*np.exp(-eta)*(eta/self.kappa - 1 - 1/self.kappa)*en/2
        numerator += csch*np.sqrt(2/self.kappa/np.pi) * \
            np.exp(-(eta**2 + self.kappa**2)/2/self.kappa)

        # Return the exact result, making use of existing relations
        gamma_reduced = self.gamma_isotensional_asymptotic_reduced(eta)
        return gamma_reduced + numerator/denominator

    def gamma_isotensional_quadrature(self, eta):
        r"""The numerical quadrature approach for the
        isotensional :math:`\gamma(\eta)`.

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Compute the nondimensional end-to-end length of the EFJC model
            at many nondimensional forces using the numerical quadrature
            approach, ensuring the exact result is always obtained:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> eta_trial = np.logspace(-2, 2, 88)
                >>> gamma_quad = model.gamma(eta_trial, approach='quadrature')
                >>> gamma_exact = model.gamma(eta_trial, approach='exact')
                >>> np.isclose(gamma_quad, gamma_exact).all()
                True

        """
        # Restrict the upper limit of integration if the link can break
        if hasattr(self, 'lambda_max'):
            upper_lim = self.lambda_max
        else:
            upper_lim = np.inf

        # Compute the rescaled partition function for a single link
        def z_fun(lambda_):
            expo_1 = eta*lambda_ - self.varepsilon*self.phi(lambda_) \
                + np.log(lambda_) - np.log(eta)
            expo_2 = expo_1 - 2*eta*lambda_
            if (np.array([expo_1, expo_2]) < self.maximum_exponent).all():
                return np.exp(expo_1) - np.exp(expo_2)
        z_rescaled = quad_vec(z_fun, self.minimum_float, upper_lim)[0]

        # Compute and return the nondimensional end-to-end length gamma
        def gamma_fun(lambda_):
            expo_1 = eta*lambda_ - self.varepsilon*self.phi(lambda_) \
                + 2*np.log(lambda_) - np.log(eta)
            expo_2 = expo_1 - 2*eta*lambda_
            expo_3 = expo_1 - np.log(lambda_) - np.log(eta)
            expo_4 = expo_3 - 2*eta*lambda_
            expos = np.array([expo_1, expo_2, expo_3, expo_4])
            if (expos < self.maximum_exponent).all():
                term_1 = np.exp(expo_1) + np.exp(expo_2)
                term_2 = np.exp(expo_3) - np.exp(expo_4)
                return (term_1 - term_2)/z_rescaled
        return quad_vec(gamma_fun, self.minimum_float, upper_lim)[0]

    def gamma_isotensional_MHMCMC(self, eta, **kwargs):
        r"""The Monte Carlo approach for the isotensional :math:`\gamma(\eta)`.

        Args:
            eta (array_like): The nondimensional force(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``MHMCMC.gamma_isotensional_MHMCMC``.

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Compute the nondimensional end-to-end length of the log-squared-FJC
            in the isotensional ensemble under a nondimensional force of 5
            using the Monte Carlo approach and ensure that it is within
            ten percent of the numerical quadrature result:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='log-squared')
                >>> gamma_MHMCMC = model.gamma(5, approach='monte carlo')
                >>> gamma_quad = model.gamma(5, approach='quadrature')
                >>> np.abs(gamma_quad - gamma_MHMCMC)/gamma_quad < 1e-1
                array([ True])

        """
        return MHMCMC(self).gamma_isotensional_MHMCMC(eta, **kwargs)

    def eta_isotensional(self, gamma, **kwargs):
        r"""Invert the isotensional :math:`\gamma(\eta)`
        for the isotensional :math:`\eta(\gamma)`.

        This function obtains the isotensional nondimensional
        single-chain mechanical response :math:`\eta(\gamma)`
        by inverting the isotensional :math:`\gamma(\eta)`.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``gamma_isotensional``.

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Check that :math:`\eta[\gamma(\eta)] = \eta\,`:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC(potential='morse')
                >>> def check_gamma(eta):
                ...     gamma_fun = lambda eta: model.gamma_isotensional(eta)
                ...     eta_fun = lambda gamma: model.eta_isotensional(gamma)
                ...     return np.isclose(eta_fun(gamma_fun(eta))[0], eta)
                >>> check_gamma(np.random.rand())
                True

        """
        def gamma_fun(eta):
            return self.gamma_isotensional(eta, **kwargs)
        return self.inv_fun_1D(gamma, gamma_fun, power=4, maxiter=250)

    def beta_Pi_config(self, config, eta_vector):
        r"""The nondimensional total potential energy of a configuration.

        This function provides the nondimensional total potential energy
        :math:`\beta \Pi = \beta U - N_b\eta\gamma`
        given the configuration of the chain, i.e. the
        vector position of each atom/hinge relative to the first one.

        Args:
            config (numpy.ndarray): The configuration of the chain,
                a :math:`(N_b+1)`-by-3 numpy array.

        Returns:
            float: The nondimensional total potential energy :math:`\beta\Pi`.

        """
        beta_U = 0
        for j in range(1, len(config)):
            lambda_ = la.norm(config[j, :] - config[j - 1, :])
            beta_U += self.beta_u(lambda_)
        return beta_U - eta_vector.dot(config[-1, :] - config[0, :])
