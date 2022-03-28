"""A module for basic single-chain model utilities.

This module consist of the class ``BasicUtility`` which contains
attributes and methods that are meant to be inherited and utilized as
basic utilities by an arbitrary single-chain model class.

Example:
    Create a freely-jointed chain single-chain model class ``FJC`` that
    inherits the ``BasicUtility`` class in order to utilize the
    Langevin function for the nondimensional isotensional mechanical response:

        >>> from ufjc.utility import BasicUtility
        >>> class FJC(BasicUtility):
        ...     def __init__(self):
        ...        BasicUtility.__init__(self)
        ...     def gamma(self, eta):
        ...         return self.langevin(eta)
        >>> model = FJC()
        >>> model.gamma([0, 0.23, 1e88])
        array([0.        , 0.07639764, 1.        ])

"""

# Import external modules
import sys
import numpy as np
from scipy.optimize import root_scalar


class BasicUtility(object):
    """The single-chain model basic utilities class.

    This class contains
    attributes and methods that are meant to be inherited and utilized as
    basic utilities by an arbitrary single-chain model class.

    Attributes:
        minimum float:
            Helps avoid overflow when dividing by small numbers.
        maximum_float:
            Helps avoid overflow in function evaluations.
        maximum_exponent:
            Helps avoid overflow in exponential functions.

    """
    def __init__(self):
        """Initializes the ``BasicUtility`` class.

        Inherits overflow-related parameters from ``sys.float_info``.
        """
        self.minimum_float = sys.float_info.min
        self.maximum_float = sys.float_info.max
        self.maximum_exponent = np.log(self.maximum_float)

    def np_array(self, input):
        """Function to return input as a numpy array.

        This function essentially serves as a wrapper for ``numpy.array``
        that returns non-array type inputs as shape (1,) numpy arrays
        rather than shape () numpy arrays, which is useful for indexing.

        Args:
            input (array_like): Anything passable to ``numpy.array``.

        Returns:
            numpy.ndarray: The input as an index-able numpy array.

        Example:
            Compare attempts to index the numpy array of an integer:

                >>> import numpy as np
                >>> from ufjc.utility import BasicUtility
                >>> try:
                ...     np.array(8)[0]
                ... except IndexError:
                ...     print('Error indexing 0-dimensional array.')
                ... finally:
                ...     BasicUtility().np_array(8)[0]
                Error indexing 0-dimensional array.
                8

        """
        if np.array(input).shape == ():
            return np.array([input])
        else:
            return np.array(input)

    def inv_fun_1D(self, y, f, power=1, guess=None, **kwargs):
        """Function to invert a mathematical function.

        This function returns the argument x given a function f(x)
        and y, the query value of the function, i.e. y = f(x).

        Note:
            This function is only meant for bijective f(x).

        Args:
            y (array_like): The value(s) of the function y.
            f (function): The function f(x).
            guess (list, optional): A list of two initial guesses.
            **kwargs: Arbitrary keyword arguments.
                Passed to ``root_scalar``.

        Returns:
            numpy.ndarray: The corresponding argument(s) x.

        Example:
            Invert a polynomial function:

                >>> from ufjc.utility import BasicUtility
                >>> f = lambda x: 1 + x**2/8 + x**3/23
                >>> BasicUtility().inv_fun_1D([55, 88], f)
                array([ 9.8712079, 11.7121826])

        """
        y = self.np_array(y)
        x = np.zeros(y.shape)
        for i in range(len(y)):

            # Function whose root will be found
            def fun(x):
                return (f(x) - y[i])**power

            # First two guesses when root-finding
            if guess is not None:
                guess_0 = guess*0.95
                guess_1 = guess
            elif y[i] == 0:
                guess_0 = 0
                guess_1 = 1e-2
            else:
                guess_0 = y[i]*0.95
                guess_1 = y[i]

            # Find the root
            x[i] = root_scalar(fun, x0=guess_0, x1=guess_1, **kwargs).root

        # Return the roots
        return x

    def log_over_sinh(self, eta):
        r"""The natural logarithm of the argument over the hyperbolic sine.

        This function avoids overflow when calculating the natural logarithm
        of the ratio of the argument to the hyperbolic sine of the argument,

        .. math::
            f(\eta) = \ln\left[\frac{\eta}{\sinh(\eta)}\right].

        It also returns zero for zero argument.

        Args:
            eta (array_like): The argument(s) of the function.

        Returns:
            numpy.ndarray: The value(s) of the function.

        Example:
            Compute the function for several argument values:

                >>> from ufjc.utility import BasicUtility
                >>> BasicUtility().log_over_sinh([0, 0.23, 888])
                array([ 0.00000000e+00, -8.80117196e-03, -8.80517881e+02])

        """
        # Determine when argument is sufficiently large
        eta = self.np_array(eta)
        eta_is_big = np.nan_to_num(eta) > 3e1
        y = np.zeros(eta.size)

        # Use asymptotic relation valid for sufficiently large arguments
        if eta_is_big.any():
            y[eta_is_big] = np.log(2*eta[eta_is_big]) - eta[eta_is_big]

        # Compute analytically otherwise, and zero where argument is zero
        where_eta_zero = eta == 0
        valid_eta = ~(eta_is_big + where_eta_zero)
        if valid_eta.any():
            y[valid_eta] = np.log(eta[valid_eta]/np.sinh(eta[valid_eta]))

        # Return the results
        return y

    def coth(self, eta):
        """The hyperbolic cotangent function.

        This function essentially serves as a wrapper for the algebraic inverse
        of ``numpy.tanh`` that avoids the singularity from evaluating at zero.

        Args:
            eta (array_like): The argument(s) of the
                hyperbolic cotangent function.

        Returns:
            numpy.ndarray: The hyperbolic cotangent(s).

        Example:
            Compute the hyperbolic cotangent for several argument values:

                >>> from ufjc.utility import BasicUtility
                >>> BasicUtility().coth([0, 0.23, np.inf])
                array([4.49423284e+307, 4.42422373e+000, 1.00000000e+000])

        """
        eta = self.np_array(eta)
        eta = np.where(eta == 0, self.minimum_float, eta)
        return 1/np.tanh(eta)

    def langevin(self, eta):
        r"""The Langevin function.

        This function computes and returns the Langevin function
        of the argument eta, which is

        .. math::
            \mathcal{L}(\eta) = \coth(\eta) - 1/\eta.

        It also avoids dividing by zero at zero argument.

        Args:
            eta (array_like): The argument(s) of the Langevin function.

        Returns:
            numpy.ndarray: The evaluated Langevin function(s).

        Example:
            Compute the Langevin function for several argument values:

                >>> from ufjc.utility import BasicUtility
                >>> BasicUtility().langevin([0, 1, 8, np.inf])
                array([0.        , 0.31303529, 0.87500023, 1.        ])

        """
        eta = self.np_array(eta)
        eta = np.where(eta == 0, self.minimum_float, eta)
        return 1/np.tanh(eta) - 1/eta

    def inv_langevin(self, gamma):
        r"""The inverse Langevin function.

        This function computes the inverse of the Langevin function
        given its query value.

        Args:
            gamma (array_like): The value(s) of the Langevin function.

        Returns:
            numpy.ndarray: The argument(s) of the Langevin function.

        Example:
            Check that :math:`\mathcal{L}^{-1}[\mathcal{L}(\eta)] = \eta\,`:

                >>> import numpy as np
                >>> from ufjc.utility import BasicUtility
                >>> def check_langevin(eta):
                ...     langevin = BasicUtility().langevin
                ...     inv_langevin = BasicUtility().inv_langevin
                ...     return np.isclose(inv_langevin(langevin(eta))[0], eta)
                >>> check_langevin(np.random.rand())
                True

        """
        return self.inv_fun_1D(gamma, self.langevin)
