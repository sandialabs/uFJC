"""The core module for the uFJC single-chain model.

This module consist of the class ``uFJC`` which, upon instantiation,
becomes a uFJC single-chain model instance with methods for computing
single-chain quantities in either thermodynamic ensemble.

Example:
    Import and create an instance of the model:

        >>> from ufjc import uFJC
        >>> model = uFJC()

"""

# Import internal modules
from .isometric import uFJCIsometric
from .potential import Potential

# Import external modules
import numpy as np
from scipy.integrate import quad


class uFJC(Potential, uFJCIsometric):
    r"""The uFJC single-chain model class.

    This class represents the uFJC single-chain model, meaning that an
    instance of this class is a uFJC single-chain model instance.
    It inherits all attributes and methods from the ``uFJCIsometric``
    class, which inherits all attributes and methods from the
    ``uFJCIsotensional`` class, which inherits all attributes and
    methods from the ``BasicUtility`` class.
    It also inherits a potential from the ``Potential`` class
    as the attribute ``pot``, a model instance itself.
    Keyword arguments are utilized during instantiation in order to
    specify model parameters; see the inherited classes and various examples
    through the documentation for more information.

    Attributes:
        N_b (int): The number of links in the chain.
        k_0 (float): The initial reaction rate coefficient.
        init_config (np.ndarray): The initial configuration for Monte Carlo.
        c_kappa (float): The constant used for ideal/Gaussian approximations.
        nondim_P_eq_normalizations (dict): The normalizations for the
            nondimensional equilibrium distribution ``nondim_P_eq``
            for each approach.
        pot (object): The link potential model instance.

    Note:
        The attributes of a instantiated model should not be changed.
        There are several quantities that are computed and stored during
        instantiation that depend on these attributes and are not re-computed.

    """
    def __init__(self, **kwargs):
        """Initializes the ``uFJC`` class,
        producing a uFJC single-chain model instance.

        First, initialize and inherit all attributes and methods
        from ``uFJCIsometric`` and ``Potential`` class instances;
        the latter assigns a potential (as attributes) to the model instance.
        Next, assign the attributes listed above.
        Then assign a default incremental link stretch function
        ``delta_lambda`` (inversion of ``eta_link``) if not already
        available from the potential choice.
        Also, adjust both ``delta_lambda`` and ``beta_u`` above
        ``eta_max`` and ``lambda_max``, if they exist for the potential,
        to prevent calculation errors.
        Finally, compute and store the normalization for ``nondim_P_eq``.

        Args:
            **kwargs: Arbitrary keyword arguments.

        """
        uFJCIsometric.__init__(self)
        Potential.__init__(self, **kwargs)

        # Get general single-chain model parameters
        self.N_b = int(kwargs.get('N_b', 1))
        self.k_0 = kwargs.get('k_0', 1)
        self.nondim_P_eq_normalizations = {}

        # Default initial configuration for Monte Carlo calculations
        self.init_config = np.append(np.array(
            [np.arange(self.N_b + 1)]).T, np.zeros((self.N_b + 1, 2)), axis=1)

        # Constant for ideal/Gaussian approximations
        self.c_kappa = \
            self.kappa*(self.kappa + 1)/(self.kappa**2 + 6*self.kappa + 3)

        # Default to inverting eta(lambda) for the link if unavailable
        if hasattr(self, 'delta_lambda') is False:
            self.delta_lambda = lambda eta: self.np_array(
               self.inv_fun_1D(eta, self.eta_link, guess=1, maxiter=250) - 1)

        # Prevent calculating bond stretch above a possible maximum force
        if hasattr(self, 'eta_max'):
            self.__old_delta_lambda = self.delta_lambda

            # Returns a high number above eta_max rather than calculating
            def delta_lambda_safe(eta):
                eta = self.np_array(eta)
                indices = eta > self.eta_max
                eta[indices] = 0
                delta_lambda = self.__old_delta_lambda(eta)
                delta_lambda[indices] = self.maximum_exponent
                return delta_lambda

            # Replace the original function with this safe function
            self.delta_lambda = delta_lambda_safe

        # Inhibit Monte Carlo sampling link stretches above lambda_max
        if hasattr(self, 'lambda_max'):
            self.__old_beta_u = self.beta_u

            # Returns a prohibitively high number above lambda_max
            def beta_u_safe(lambda_):
                return self.__old_beta_u(lambda_) + \
                    self.maximum_exponent * \
                    np.heaviside(lambda_ - self.lambda_max, 0)

                # Replace the original function with this safe function
            self.beta_u = beta_u_safe

    def gamma(self, eta, **kwargs):
        r"""The nondimensional end-to-end length.

        This function computes the scalar nondimensional end-to-end length as a
        function of the scalar nondimensional force, :math:`\gamma(\eta)`.
        In the isotensional ensemble, it is given by
        :cite:`fiasconaro2019analytical`

        .. math::
            \gamma(\eta) = \frac{\partial}{\partial\eta}\,\ln\mathfrak{z}(\eta)
            ,

        where the single-link (since :math:`\gamma(\eta)` is independent of
        :math:`N_b` in the isotensional ensemble for the uFJC model)
        nondimensional isotensional partition function is given by
        :cite:`buche2022on`

        .. math::
            \mathfrak{z}(\eta) =
            4\pi\int \frac{\sinh(\lambda\eta)}{\lambda\eta} \,
            e^{-\varepsilon\phi(\lambda)} \, \lambda^2 d\lambda
            .

        This function is rarely exactly
        analytically available, is sometimes accurately approximated using
        analytically tractable asymptotic relations, and is otherwise
        numerically calculated using Monte Carlo or quadrature approaches.
        For the uFJC model, accurate asymptotic relations are available,
        and for the EFJC model (harmonic link potentials), it is actually
        possible to exactly find :math:`\gamma(\eta)` closed-form
        :cite:`balabaev2009extension`.

        Args:
            eta (array_like): The nondimensional force(s).
            **kwargs: Arbitrary keyword arguments.

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Examples:
            Create an EFJC model with a nondimensional link energy
            :math:`\varepsilon=23`, and calculate the nondimensional
            end-to-end length under an nondimensional force of 8
            using many of the available approaches:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='harmonic', varepsilon=23)
                >>> model.gamma(8, approach='exact')
                array([1.25508427])
                >>> model.gamma(8, approach='asymptotic')
                array([1.25508427])
                >>> model.gamma(8, approach='reduced')
                array([1.22282631])
                >>> model.gamma(8, approach='quadrature')
                array([1.25508427])

        Example:
            Verify that the ideal approximation is valid
            for small nondimensional forces:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> gamma_exact = model.gamma(1e-2, approach='exact')
                >>> gamma_ideal = model.gamma(1e-2, ideal=True)
                >>> np.isclose(gamma_exact, gamma_ideal)
                array([ True])

        Example:
            Verify that the single-chain mechanical response is the same
            in the isotensional ensemble as it is in the isometric ensemble
            when utilizing the Legendre transformation method:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='harmonic', varepsilon=23)
                >>> eta = np.random.rand(88)
                >>> (np.abs(model.gamma(eta, ensemble='isotensional') -
                ...         model.gamma(eta, ensemble='isometric',
                ...                     method='legendre')) < 1e-6).all()
                True

        """
        if kwargs.get('ideal', False) is True:
            return self.np_array(eta)/(3*self.c_kappa)
        else:
            if kwargs.get('ensemble', 'isotensional') == 'isometric':
                return self.gamma_isometric(eta, **kwargs)
            else:
                return self.gamma_isotensional(eta, **kwargs)

    def eta(self, gamma, **kwargs):
        r"""The nondimensional force.

        This function computes the scalar nondimensional force as a function of
        the scalar nondimensional end-to-end length, :math:`\eta(\gamma)`.
        In the isometric ensemble, it is given by :cite:`manca2012elasticity`

        .. math::
            \gamma(\eta) = -\frac{1}{N_b}
            \frac{\partial}{\partial\gamma}
            \,\ln\mathfrak{Q}(\boldsymbol{\gamma})
            ,

        where the nondimensional isometric partition function is given by
        :cite:`buche2021fundamental`

        .. math::
            \mathfrak{Q}(\boldsymbol{\gamma}) =
            \iiint d^3\boldsymbol{\lambda}_1
            \,e^{-\varepsilon\phi(\lambda_1)} \cdots
            \iiint d^3\boldsymbol{\lambda}_{N_b}
            \,e^{-\varepsilon\phi(\lambda_{N_b})}
            ~\delta^3\left(\sum_{j=1}^{N_b}\boldsymbol{\lambda}_j
            - \boldsymbol{\gamma}\right)
            .

        This is analytically intractable for nearly all single-chain models,
        including the uFJC model. Fortunately, there are both exact numerical
        methods, such as Monte Carlo calculations, and accurate approximation
        methods, such as the Legendre transformation method which becomes
        valid relatively quickly as :math:`N_b` increases.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Compute the nondimensional force at an end-to-end length in the
            isometric ensemble, which (by default) is done using the Legendre
            transformation method and the asymptotic approach to get
            :math:`\gamma(\eta)` in the isotensional ensemble. Using the
            Legendre transformation method with a given isotensional approach
            for the isometric ensemble is equivalent to using the same approach
            in the isotensional ensemble:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='harmonic', varepsilon=23, N_b=8)
                >>> model.eta(1.25508427, ensemble='isometric')
                array([8.00000001])
                >>> model.eta(1.25508427, ensemble='isotensional')
                array([8.00000001])

        Example:
            Verify that the ideal approximation is valid for small
            nondimensional end-to-end lengths:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> eta_exact = model.eta(1e-2, approach='exact')
                >>> eta_ideal = model.eta(1e-2, ideal=True)
                >>> np.abs((eta_exact - eta_ideal)/eta_exact) < 1e-3
                array([ True])

        """
        if kwargs.get('ideal', False) is True:
            return 3*self.c_kappa*self.np_array(gamma)
        else:
            if kwargs.get('ensemble', 'isotensional') == 'isometric':
                return self.eta_isometric(gamma, **kwargs)
            else:
                return self.eta_isotensional(gamma, **kwargs)

    def vartheta(self, gamma, **kwargs):
        r"""The nondimensional Helmholtz free energy per link.

        This function compute the nondimensional Helmholtz free energy per link
        in the isometric ensemble, which is generally given by

        .. math::
            \vartheta(\boldsymbol{\gamma}) =
            -\ln\left[\mathfrak{Q}(\boldsymbol{\gamma})\right]^{1/N_b}
            .

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``vartheta_isometric``.

        Returns:
            numpy.ndarray: The nondimensional Helmholtz
            free energy(s) per link.

        Example:
            Create a uFJC single-chain model with log-squared link potentials
            and calculate the Helmholtz free energy per link at
            a nondimensional end-to-end length of a half:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='log-squared')
                >>> model.vartheta(0.5, method='legendre', approach='reduced')
                array([0.39097337])

        Example:
            Verify that the ideal approximation is valid for small
            nondimensional end-to-end lengths:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> vartheta = model.vartheta(1e-2, approach='reduced')
                >>> vartheta_ideal = model.vartheta(1e-2, ideal=True)
                >>> np.abs((vartheta - vartheta_ideal)/vartheta) < 1e-1
                array([ True])

        """
        if kwargs.get('ideal', False) is True:
            return 1.5*self.c_kappa*self.np_array(gamma)**2
        else:
            return self.vartheta_isometric(gamma, **kwargs)

    def nondim_P_eq(self, gamma, **kwargs):
        r"""The nondimensional equilibrium distribution function.

        This function computes the nondimensional probability density
        distribution of nondimensional end-to-end vectors at equilibrium,

        .. math::
            \mathscr{P}_\mathrm{eq}(\boldsymbol{\gamma}) =
            \frac{e^{-N_b\,\vartheta(\boldsymbol{\gamma})}}
            {\iiint e^{-N_b\,\vartheta(\tilde{\boldsymbol{\gamma}})}
            \,d^3\tilde{\boldsymbol{\gamma}}}
            .

        Since the uFJC model is spherically-symmetric in
        :math:`\boldsymbol{\gamma}`,  ``nondim_P_eq`` is only a function
        of the scalar nondimensional end-to-end length :math:`\gamma`,
        and the normalization can be calculated
        using a one-dimensional integral over :math:`\gamma` instead
        :cite:`buche2020statistical`,

        .. math::
            \iiint e^{-N_b\,\vartheta(\tilde{\boldsymbol{\gamma}})}
            \,d^3\tilde{\boldsymbol{\gamma}} =
            4\pi\int e^{-N_b\,\vartheta(\tilde{\gamma})}
            \,\tilde{\gamma}^2 d\tilde{\gamma}
            .

        The normalization for the distribution is computed for a given approach
        once the function is called using that approach.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments. Passed to ``vartheta``.

        Returns:
            numpy.ndarray: The nondimensional equilibrium
            probability density(s).

        Example:
            Evaluate the probability density at an end-to-end length of 0.23
            using two different approximation methods:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='morse', N_b=8)
                >>> model.nondim_P_eq(23e-2)
                array([4.19686303])
                >>> model.nondim_P_eq(23e-2, gaussian=True)
                array([3.86142625])

        Example:
            Verify that the probability density is normalized:

                >>> from ufjc import uFJC
                >>> model = uFJC()
                >>> from scipy.integrate import quad
                >>> integrand = lambda gamma: \
                ...     4*np.pi*gamma**2*model.nondim_P_eq(gamma)
                >>> P_tot_eq = quad(integrand, 0, np.inf)[0]
                >>> np.isclose(P_tot_eq, 1)
                True

        """
        if kwargs.get('gaussian', False) is True:
            a = 1.5*self.N_b*self.c_kappa
            return (a/np.pi)**(1.5)*np.exp(-a*self.np_array(gamma)**2)
        else:

            # A variable string based on the approach
            var_str = 'nondim_P_eq_normalization_' + \
                str(kwargs.get('approach', 'asymptotic'))

            # Compute the normalization if not done already for the approach
            if self.nondim_P_eq_normalizations.get(var_str) is None:

                # Set upper limit of integration
                if hasattr(self, 'lambda_max') is True:
                    upper_lim = self.lambda_max
                else:
                    upper_lim = np.inf

                # Compute the normalzation and attribute it
                self.nondim_P_eq_normalizations[var_str] = \
                    quad(lambda gamma: 4*np.pi*gamma**2 *
                         np.exp(-self.N_b*self.vartheta(gamma)),
                         self.minimum_float, upper_lim)[0]

            # Return the nondimensional equilibrium probability density
            return np.exp(-self.N_b*self.vartheta(gamma, **kwargs)) / \
                self.nondim_P_eq_normalizations[var_str]

    def nondim_g_eq(self, gamma, **kwargs):
        r"""The nondimensional equilibrium radial distribution function.

        This function computes the nondimensional radial probability density
        distribution of nondimensional end-to-end lengths at equilibrium,

        .. math::
            \mathscr{g}_\mathrm{eq}(\gamma) =
            4\pi\gamma^2 \mathscr{P}_\mathrm{eq}(\gamma)
            .

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).
            **kwargs: Arbitrary keyword arguments. Passed to ``nondim_P_eq``.

        Returns:
            numpy.ndarray: The nondimensional equilibrium radial
            probability density(s).

        Example:
            Verify that the function is normalized by integrating over all
            permissible end-to-end lengths, which in effect ensures that the
            total probability is unity:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> model = uFJC(potential='morse')
                >>> from scipy.integrate import quad
                >>> P_tot_eq = quad(model.nondim_g_eq, 0, model.lambda_max)[0]
                >>> np.isclose(P_tot_eq, 1)
                True

        """
        return 4*np.pi*self.np_array(gamma)**2 * \
            self.nondim_P_eq(gamma, **kwargs)

    def k(self, gamma):
        r"""The net forward reaction rate coefficient function.

        This function computes the net forward reaction rate coefficient
        as a function of the nondimensional end-to-end length,
        i.e. :math:`k(\gamma)=N_bk'(\gamma)`, where :math:`k'(\gamma)` is the
        forward reaction rate coefficient function for a single link.
        This function is obtained using the axioms of transition state theory
        :cite:`zwanzig2001nonequilibrium`,
        and is currently implemented to use the Legendre transformation method
        and the reduced asymptotic approach
        :cite:`buche2021chain`.

        Args:
            gamma (array_like): The nondimensional end-to-end length(s).

        Returns:
            numpy.ndarray: The net forward reaction rate coefficient(s).

        Example:
            Compute the net forward reaction rate coefficients at several
            nondimensional end-to-end lengths:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='morse')
                >>> model.k([0.8, 1.01, 1.6])
                array([6.23755528e-01, 9.80472713e-01, 6.59105919e+07])

        """
        # Compute the corresponding force and link stretch
        eta = self.np_array(self.eta(gamma))
        lambda_ = 1 + self.delta_lambda(eta)

        # Compute the free energy barrier for a link, translated so k(0)=k_0
        beta_delta_Psi_TS_single_link = \
            self.varepsilon*(self.phi(1) - self.phi(lambda_))
        if hasattr(self, 'lambda_max'):
            beta_delta_Psi_TS_single_link += -self.log_over_sinh(eta) + \
                self.log_over_sinh(self.lambda_max*eta) + \
                self.lambda_max*eta*self.coth(self.lambda_max*eta) + \
                - eta*self.coth(eta)

        # Avoid overflow when returning the results
        exponent = np.log(self.N_b*self.k_0) - beta_delta_Psi_TS_single_link
        exponent[exponent > self.maximum_exponent] = self.maximum_exponent
        return np.exp(exponent)
