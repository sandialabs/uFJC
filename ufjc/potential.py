"""A module for handling potentials.

This module contains several different classes representing potentials,
each having methods to compute relevant nondimensional quantities as
functions of nondimensional force or stretch.
This module also contains the parent class ``Potential`` that is used to
assign a potential to a given model using keyword arguments.

Examples:
    Create a Lennard-Jones potential model with
    a nondimensional potential energy scale of 8 and evaluate
    the nondimensional potential energy at a stretch of 1.23:

        >>> from ufjc.potential import LennardJonesPotential
        >>> model = LennardJonesPotential(varepsilon=8)
        >>> model.beta_u(1.23)
        -3.953345685631384

    Create a single-link model in one dimension, instantiate it with
    the Morse potential, and compute the incremental link stretch under
    a nondimensional force of 8:

        >>> from ufjc.potential import Potential
        >>> class Link1D(Potential):
        ...     def __init__(self, **kwargs):
        ...         Potential.__init__(self, **kwargs)
        >>> Link1D(potential='morse').delta_lambda(8)
        0.04890980361596759

"""

# Import external modules
import numpy as np
from scipy.special import lambertw


class HarmonicPotential(object):
    r"""The harmonic potential.

    Attributes:
        varepsilon (float): The nondimensional energy scale.
        kappa (float): The nondimensional stiffness :math:`\kappa=\varepsilon`.
        c (float): The correction parameter :math:`c=1`.

    """
    def __init__(self, **kwargs):
        """Initializes the ``HarmonicPotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Can be used to specify ``varepsilon`` (default 88).

        """
        self.varepsilon = kwargs.get('varepsilon', 88)
        self.kappa = self.varepsilon
        self.c = 1

    def phi(self, lambda_):
        r"""The scaled nondimensional potential energy function,

            .. math::
                \phi(\lambda) = \frac{1}{2}(\lambda-1)^2
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The scaled nondimensional potential energy(s).

        """
        return 0.5*(lambda_ - 1)**2

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)

    def eta_link(self, lambda_):
        r"""The nondimensional force as a function of stretch,

            .. math::
                \eta(\lambda) = \varepsilon(\lambda - 1)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Compute the nondimensional force at a sample stretch:

                >>> from ufjc.potential import HarmonicPotential
                >>> HarmonicPotential().eta_link(1.8)
                70.4

        """
        return self.varepsilon*(lambda_ - 1)

    def delta_lambda(self, eta):
        r"""The incremental stretch as a function of nondimensional force,

            .. math::
                \Delta\lambda(\eta) = \frac{\eta}{\varepsilon}
                .

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The incremental stretch(s).

        """
        return eta/self.varepsilon


class LogSquaredPotential(object):
    r"""The log-squared potential :cite:`mao2017rupture`.

    Attributes:
        varepsilon (float): The nondimensional energy scale.
        kappa (float): The nondimensional stiffness :math:`\kappa=\varepsilon`.
        c (float): The correction parameter :math:`c=2/5`.
        eta_max (float): The maximum nondimensional force
            :math:`\eta_\mathrm{max} = e^{-1}\varepsilon`.
        lambda_max (float): The stretch at the maximum nondimensional force,
            :math:`\lambda_\mathrm{max} = e^{1}`.

    """
    def __init__(self, **kwargs):
        """Initializes the ``LogSquaredPotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Can be used to specify ``varepsilon`` (default 88).

        """
        self.varepsilon = kwargs.get('varepsilon', 88)
        self.kappa = self.varepsilon
        self.c = 2/5
        self.eta_max = self.varepsilon/np.exp(1)
        self.lambda_max = np.exp(1)

    def phi(self, lambda_):
        r"""The scaled nondimensional potential energy function,

            .. math::
                \phi(\lambda) = \frac{1}{2}\big[\ln(\lambda)\big]^2
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The scaled nondimensional potential energy(s).

        """
        return 0.5*np.log(lambda_)**2

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)

    def eta_link(self, lambda_):
        r"""The nondimensional force as a function of stretch,

            .. math::
                \eta(\lambda) = \varepsilon\,\frac{\ln(\lambda)}{\lambda}
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Compute the nondimensional force at a sample stretch:

                >>> from ufjc.potential import LogSquaredPotential
                >>> LogSquaredPotential().eta_link(1.8)
                28.736236950770266

        """
        return self.varepsilon*np.log(lambda_)/lambda_

    def delta_lambda(self, eta):
        r"""The incremental stretch as a function of nondimensional force,

            .. math::
                \Delta\lambda(\eta) = e^{-\mathcal{W}(-\eta/\varepsilon)}
                ,\qquad \eta\in[0,\eta_\mathrm{max}]
                .

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The incremental stretch(s).

        """
        return (np.exp(-lambertw(-eta/self.varepsilon)) - 1).real


class MorsePotential(object):
    r"""The Morse potential :cite:`morse1929diatomic`.

    Attributes:
        varepsilon (float): The nondimensional energy scale.
        alpha (float): The Morse parameter.
        kappa (float): The nondimensional stiffness
            :math:`\kappa=2\varepsilon\alpha^2`.
        c (float): The correction parameter
            :math:`c=1/(1+3\alpha/2)`.
        eta_max (float): The maximum nondimensional force
            :math:`\eta_\mathrm{max} = \sqrt{\kappa\varepsilon/8}`.
        lambda_max (float): The stretch at the maximum nondimensional force,
            :math:`\lambda_\mathrm{max} = 1+\ln(2)/\alpha`.

    """
    def __init__(self, **kwargs):
        """Initializes the ``MorsePotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Can be used to specify ``varepsilon`` (default 88)
                and ``alpha`` (default 1).

        """
        self.varepsilon = kwargs.get('varepsilon', 88)
        self.alpha = kwargs.get('alpha', 1)
        self.kappa = 2*self.varepsilon*self.alpha**2
        self.c = 1/(1 + 3/2*self.alpha)
        self.eta_max = np.sqrt(self.kappa*self.varepsilon/8)
        self.lambda_max = 1 + np.log(2)/self.alpha

    def phi(self, lambda_):
        r"""The scaled nondimensional potential energy function,

            .. math::
                \phi(\lambda) = \left[1
                    - e^{-\alpha(\lambda - 1)}\right]^2
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The scaled nondimensional potential energy(s).

        """
        return (1 - np.exp(-self.alpha*(lambda_ - 1)))**2

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)

    def eta_link(self, lambda_):
        r"""The nondimensional force as a function of stretch,

            .. math::
                \eta(\lambda) = 2\alpha\varepsilon e^{-\alpha(\lambda - 1)}
                    \left[1 - e^{-\alpha(\lambda - 1)}\right]^2
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional force(s).

        Example:
            Compute the nondimensional force at a sample stretch:

                >>> from ufjc.potential import MorsePotential
                >>> MorsePotential().eta_link(1.23)
                28.731992431367807

        """
        return 2*self.alpha*self.varepsilon * \
            np.exp(-self.alpha*(lambda_ - 1)) * \
            (1 - np.exp(-self.alpha*(lambda_ - 1)))

    def delta_lambda(self, eta):
        r"""The incremental stretch as a function of nondimensional force,

            .. math::
                \Delta\lambda(\eta) =
                    \ln\left(\frac{2}{1 + \sqrt{1 -
                        \eta/\eta_\mathrm{max}}}\right)^{1/\alpha}
                ,\qquad \eta\in[0,\eta_\mathrm{max}]
                .

        Args:
            eta (array_like): The nondimensional force(s).

        Returns:
            numpy.ndarray: The incremental stretch(s).

        """
        return np.log(2/(1 + np.sqrt(1 - eta/self.eta_max)))/self.alpha


class LennardJonesPotential(object):
    r"""The Lennard-Jones potential :cite:`jones1924determinationii`.

    Attributes:
        varepsilon (float): The nondimensional energy scale.
        kappa (float): The nondimensional stiffness
            :math:`\kappa=72\varepsilon`.
        c (float): The correction parameter
            :math:`c=2/23`.
        eta_max (float): The maximum nondimensional force
            :math:`\eta_\mathrm{max} = \eta(\lambda_\mathrm{max})`.
        lambda_max (float): The stretch at the maximum nondimensional force,
            :math:`\lambda_\mathrm{max} = (13/7)^{1/6}`.

    """
    def __init__(self, **kwargs):
        """Initializes the ``LennardJonesPotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Can be used to specify ``varepsilon`` (default 88).

        """
        self.varepsilon = kwargs.get('varepsilon', 88)
        self.kappa = 72*self.varepsilon
        self.c = 2/23
        self.lambda_max = (13/7)**(1/6)
        self.eta_max = self.eta_link(self.lambda_max)

    def phi(self, lambda_):
        r"""The scaled nondimensional potential energy function,

            .. math::
                \phi(\lambda) =
                    \frac{1}{\lambda^{12}} - \frac{2}{\lambda^6}
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The scaled nondimensional potential energy(s).

        """
        return 1/lambda_**12 - 2/lambda_**6

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)

    def eta_link(self, lambda_):
        r"""The nondimensional force as a function of stretch,

            .. math::
                \eta(\lambda) = 12\varepsilon\left(
                    \frac{1}{\lambda^7} - \frac{1}{\lambda^{13}}\right)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional force(s).

        """
        return self.varepsilon*(12/lambda_**7 - 12/lambda_**13)


class MiePotential(object):
    r"""The Mie potential :cite:`mie1903kinetischen`.

    Attributes:
        varepsilon (float): The nondimensional energy scale.
        n (float): The repulsive exponent.
        m (float): The attractive exponent.
        kappa (float): The nondimensional stiffness
            :math:`\kappa=nm\varepsilon`.
        c (float): The correction parameter
            :math:`c=\frac{4m(m+1)-2n(n+1)}{2m(m^2+5m+4)-n(n^2+5n+4)}`.
        eta_max (float): The maximum nondimensional force
            :math:`\eta_\mathrm{max} = \eta(\lambda_\mathrm{max})`.
        lambda_max (float): The stretch at the maximum nondimensional force,
            :math:`\lambda_\mathrm{max}=\left[\frac{n+1}{m+1}\right]^{1/(n-m)}`
            .

    """
    def __init__(self, **kwargs):
        """Initializes the ``MiePotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Can be used to specify ``varepsilon`` (default 88)
                ``n`` (default 12) and ``m`` (default 6).

        """
        self.varepsilon = kwargs.get('varepsilon', 88)
        self.m = kwargs.get('m', 6)
        self.n = kwargs.get('n', 12)
        self.kappa = self.n*self.m*self.varepsilon
        self.c = (4*self.m*(self.m + 1) - 2*self.n*(self.n + 1))/(
            2*self.m*(self.m**2 + 5*self.m + 4) -
            self.n*(self.n**2 + 5*self.n + 4))
        self.lambda_max = ((self.n + 1)/(self.m + 1))**(1/(self.n - self.m))
        self.eta_max = self.eta_link(self.lambda_max)

    def phi(self, lambda_):
        r"""The scaled nondimensional potential energy function,

            .. math::
                \phi(\lambda) = \frac{1}{(n-m)}
                    \left(\frac{m}{\lambda^n} - \frac{n}{\lambda^m}\right)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The scaled nondimensional potential energy(s).

        """
        return (self.m/lambda_**self.n -
                self.n/lambda_**self.m)/(self.n - self.m)

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)

    def eta_link(self, lambda_):
        r"""The nondimensional force as a function of stretch,

            .. math::
                \eta(\lambda) = \varepsilon\, \frac{nm}{(n-m)}
                    \left(\frac{1}{\lambda^{m+1}}
                        - \frac{1}{\lambda^{n+1}}\right)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional force(s).

        """
        return self.varepsilon*self.n*self.m/(self.n - self.m)*(
            1/lambda_**(self.m + 1) - 1/lambda_**(self.n + 1))


class PolynomialPotential(object):
    r"""A polynomial potential.

    Attributes:
        varepsilon (float): The nondimensional energy scale.
        coefficients (array_like): The coefficients :math:`a_k`.
        kappa (float): The nondimensional stiffness
            :math:`\kappa=a_1\varepsilon`.
        c (float): The correction parameter
            :math:`c=(1 - \frac{a_2}{2a_1})^{-1}`.

    """
    def __init__(self, **kwargs):
        """Initializes the ``PolynomialPotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Can be used to specify ``varepsilon`` (default 88)
                and ``coefficients`` (default 1).

        """
        self.varepsilon = kwargs.get('varepsilon', 88)
        coef = np.array(kwargs.get('coefficients', [1, 0]))
        self.eta_c = np.append(np.array([0]), coef)
        self.phi_c = np.append([0], self.eta_c) * \
            np.append([0, 0], [1/n for n in range(2, len(coef) + 2)])
        self.kappa = self.varepsilon*self.eta_c[1]
        if len(self.phi_c) > 3:
            self.c = 1/(1 - self.phi_c[3]/self.phi_c[2]/2)

    def phi(self, lambda_):
        r"""The scaled nondimensional potential energy function,

            .. math::
                \phi(\lambda) = \sum_{k=2} \frac{a_{k-1}}{k}\, \lambda^k
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The scaled nondimensional potential energy(s).

        """
        return np.polynomial.Polynomial(self.phi_c)(lambda_ - 1)

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)

    def eta_link(self, lambda_):
        r"""The nondimensional force as a function of stretch,

            .. math::
                \eta(\lambda) = \varepsilon\sum_{k=1} a_k\, \lambda^k
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional force(s).

        """
        return self.varepsilon * \
            np.polynomial.Polynomial(self.eta_c)(lambda_ - 1)


class CustomPotential(object):
    r"""A custom user-defined potential.

    Attributes:
        potential (str): The potential name.
        varepsilon (float): The nondimensional energy scale.
        phi (function): The scaled nondimensional energy function.
        eta_link (function): The nondimensional force as a function of stretch.
        delta_lambda (function): The incremental stretch as a function of the
            nondimensional force (optional).
        kappa (float): The nondimensional stiffness.
        c (float): The correction parameter.

    """
    def __init__(self, **kwargs):
        """Initializes the ``CustomPotential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.

        """
        self.potential = kwargs.get('potential', 'custom')
        self.varepsilon = kwargs.get('varepsilon')
        self.phi = kwargs.get('phi')
        self.eta_link = kwargs.get('eta_link')
        self.kappa = kwargs.get('kappa')
        self.c = kwargs.get('c')

        # No defaults
        if 'delta_lambda' in kwargs:
            self.delta_lambda = kwargs.get('delta_lambda')
        if 'lambda_max' in kwargs:
            self.lambda_max = kwargs.get('lambda_max')
        if 'eta_max' in kwargs:
            self.eta_max = kwargs.get('eta_max')

    def beta_u(self, lambda_):
        r"""The nondimensional potential energy function,

            .. math::
                \beta u(\lambda) = \varepsilon\phi(\lambda)
                .

        Args:
            lambda_ (array_like): The stretch(s).

        Returns:
            numpy.ndarray: The nondimensional potential energy(s).

        """
        return self.varepsilon*self.phi(lambda_)


class Potential(object):
    r"""A class to assign a potential to a given model through inheritance.

    Attributes:
        potential (str): The potential type.
        pot (object): The potential model instance.

    """
    def __init__(self, **kwargs):
        """Initializes the ``Potential`` class.

        Args:
            **kwargs: Arbitrary keyword arguments.
                Used to specify the potential and then passed to
                an instantiation of that potential.

        Note:
            An improperly-specified potential will default to harmonic:

                >>> from ufjc import uFJC
                >>> model = uFJC(potential='blah')
                Potential "blah" invalid, defaulting to "harmonic" potential.

        """
        self.potential = kwargs.get('potential', 'harmonic')
        if self.potential == 'harmonic':
            self.pot = HarmonicPotential(**kwargs)
        elif self.potential == 'log-squared':
            self.pot = LogSquaredPotential(**kwargs)
        elif self.potential == 'morse':
            self.pot = MorsePotential(**kwargs)
        elif self.potential == 'lennard-jones':
            self.pot = LennardJonesPotential(**kwargs)
        elif self.potential == 'mie':
            self.pot = MiePotential(**kwargs)
        elif self.potential == 'polynomial':
            self.pot = PolynomialPotential(**kwargs)
        elif self.potential == 'custom':
            self.pot = CustomPotential(**kwargs)
        else:
            print('Potential "' + self.potential +
                  '" invalid, defaulting to "harmonic" potential.')
            self.pot = HarmonicPotential(**kwargs)

    def __getattr__(self, attr):
        """Inherit attributes and methods from the chosen potential.

            Note:
                When accessing an attribute of the model to which the potential
                is assigned, but no attribute is found, and ``AttributeError``
                will be raised with respect to the ``pot`` instance, not the
                overall model instance. This is sub-optimal but typically
                acceptable, but still should be repaired in the future.

        """
        return getattr(self.pot, attr)
