"""A module for the Metropolis-Hastings Markov chain Monte Carlo method.

This module consist of the class ``MHMCMC`` which contains
attributes and methods that allow single-chain quantities
to be computed using the Metropolis-Hastings Markov chain Monte Carlo method.

Example:
    Import and instantiate the class with a single-chain model:

        >>> from ufjc import uFJC
        >>> from ufjc.monte_carlo import MHMCMC
        >>> class_instance = MHMCMC(uFJC())

"""

# Import internal modules
from .utility import BasicUtility

# Import external modules
import numpy as np
import numpy.linalg as la
import multiprocessing as mp


class MHMCMC(BasicUtility):
    """The Metropolis-Hastings Markov chain Monte Carlo class.

    This class contains attributes and methods used to compute
    statistical thermodynamic quantities using the
    Metropolis-Hastings Markov chain Monte Carlo method
    :cite:`haile1992molecular` for a given single-chain model.

    Attributes:
        N_b:
            The number of chain links.
            Inherited from ``single_chain_model``.
        beta_U_config:
            The configurational potential energy.
            Inherited from ``single_chain_model``.
        beta_Pi_config:
            The configurational total potential energy.
            Inherited from ``single_chain_model``.
        init_config:
            The default initial configuration.
            Inherited from ``single_chain_model``.

    """
    def __init__(self, single_chain_model):
        """Initializes the ``MHMCMC`` class.

        Initialize and inherit all attributes and methods
        from a ``BasicUtility`` class instance, and certain
        attributes and methods from the given single-chain model.

        Args:
            single_chain_model(object): The single-chain model instance.

        """
        BasicUtility.__init__(self)

        # Inherit necessary attributes from the single-chain model
        self.N_b = single_chain_model.N_b
        self.init_config = single_chain_model.init_config
        self.beta_U_config = single_chain_model.beta_U_config
        self.beta_Pi_config = single_chain_model.beta_Pi_config

    def trial_config(self, prev_config, **kwargs):
        """Generates trial configurations given a previous configuration.

        This function returns a trial configuration given a previous starting
        configuration.
        The possibility of rigid-body translation is prevented
        by fixing the first coordinate at the origin.
        In the isometric ensemble, the last coordinate is left unchanged
        in order to keep the end-to-end vector fixed.

        Args:
            prev_config (numpy.ndarray): The previous configuration.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            numpy.ndarray: The generated trial configuration.

        Example:
            Generate a trial configuration in the isometric ensemble and
            ensure that the end-to-end remains unchanged:

                >>> from ufjc import uFJC
                >>> from ufjc.monte_carlo import MHMCMC
                >>> model = uFJC(N_b=8)
                >>> prev_config = model.init_config
                >>> trial_config = MHMCMC(model).trial_config(prev_config,
                ...     ensemble='isometric')
                >>> (prev_config == trial_config)[(0,-1),:]
                array([[ True,  True,  True],
                       [ True,  True,  True]])

        """
        cov_config = kwargs.get('transition_variance', 1e-2)*np.eye(3)
        delta_config = np.random.multivariate_normal(
            np.zeros(3), cov_config, len(prev_config))
        delta_config[0, :] = 0
        if kwargs.get('ensemble', 'isotensional') == 'isometric':
            delta_config[-1, :] = 0
        return prev_config + delta_config

    def mh_next_config(self, prev_config, **kwargs):
        """Generates accepted configurations given a previous configuration.

        This function returns an accepted configuration given a
        previous configuration. Trial configurations are generated and
        repeatedly rejected until accepted as the next configuration
        based on the Metropolis-Hastings algorithm.

        Args:
            prev_config (numpy.ndarray): The previous configuration.
            **kwargs: Arbitrary keyword arguments.
                Passed to ``trial_config``.

        Returns:
            numpy.ndarray: The accepted (next) configuration.

        """
        next_config = None
        urng = kwargs.get('urng', np.random.uniform)
        while next_config is None:
            trial_config = self.trial_config(prev_config, **kwargs)
            delta_bEc = self.__bEc(trial_config) - self.__bEc(prev_config)
            if (delta_bEc < 0) or (np.exp(-delta_bEc) > urng()):
                next_config = trial_config
        return next_config

    def mh_next_config_append(self, prev_configs, **kwargs):
        """Append next configuration onto a history of previous configurations.

        This function calls ``mh_next_config`` and appends the resulting next
        configuration onto a given history of previous configurations.

        Args:
            prev_configs (numpy.ndarray): The previous
                history of configurations.
            **kwargs: Arbitrary keyword arguments.
                Passed to ``mh_next_config``.

        Returns:
            numpy.ndarray: The updated history of configurations.

        """
        return np.append(prev_configs, np.array(
            [self.mh_next_config(prev_configs[-1], **kwargs)]), axis=0)

    def burn_in(self, init_config, **kwargs):
        """Burn-in (throw away samples) until in the effective sampling region.

        This function runs a Monte Carlo calculation without utilizing samples,
        instead throwing them away until a desired number have been thrown away
        or a desired tolerance has been reached.
        This is typically done to obtain a burned-in configuration to use as
        the initial configuration for an actual Monte Carlo calculation.

        Args:
            init_config (numpy.ndarray):
                The initial configuration.
            **kwargs: Arbitrary keyword arguments.
                Passed to ``mh_next_config_append``.

        Returns:
            numpy.ndarray: The final, burned-in configuration.

        """
        config_history = np.array([init_config])
        prev_rmc = init_config

        # Generate new configurations until certain criterion are met
        tol = kwargs.get('tol', np.inf)
        num_samples_burn_in = int(kwargs.get('num_samples_burn_in', int(1e3)))
        n_res = tol + 1
        while (n_res > tol) or (len(config_history) < num_samples_burn_in):

            # Generate an accepted configuration
            config_history = self.mh_next_config_append(
                config_history, **kwargs)

            # Compute the normalized residual in mean configuration
            next_rmc = np.mean(config_history, axis=0)
            n_res = la.norm(next_rmc - prev_rmc)/la.norm(prev_rmc)
            prev_rmc = next_rmc

        # Return the configuration at the end of the burn-in period
        return config_history[-1]

    def parallel_calculation(self, serial_fun, **kwargs):
        """Monte Carlo calculation averaged over several parallel processes.

        This function performs several Monte Carlo calculations in parallel
        and returns the average result from all processes.
        Each Monte Carlo calculation is performed using ``serial_fun``.
        The default is to utilize all available processors.

        Args:
            serial_fun(function):
                The function for a single Monte Carlo calculation.
            **kwargs: Arbitrary keyword arguments.
                Passed to ``serial_fun``.

        Returns:
            float: The process-averaged result.

        """
        burned_in_config = self.burn_in(self.init_config, **kwargs)
        num_processes = kwargs.get('num_processes', mp.cpu_count())
        if num_processes > 1:

            # Create an output queue
            output = mp.Queue()

            # Set up the function run by each process
            def fun(seed, output):
                output.put(serial_fun(burned_in_config,
                           urng=np.random.RandomState(seed).random, **kwargs))

            # Create the processes, each with a random seed
            processes = [mp.Process(target=fun, args=(seed, output))
                         for seed in np.random.randint(88, size=num_processes)]

            # Run each process
            for p in processes:
                p.start()

            # Exit each process after running them all
            for p in processes:
                p.join()

            # Get the results from each process from the output queue
            process_results = [output.get() for p in processes]

            # Calculate and return the mean over all processes
            return np.mean(process_results)

        # Otherwise do a single serial calculation
        else:
            return serial_fun(burned_in_config, **kwargs)

    def gamma_isotensional_MHMCMC_serial(self, config, **kwargs):
        r"""Serial Monte Carlo calculation function for
        the isotensional :math:`\gamma(\eta)`.

        This function computes the nondimensional end-to-end length
        using a Monte Carlo calculation given an initial configuration.

        Args:
            config (numpy.ndarray):
                The initial configuration, typically already burned-in.
            **kwargs: Arbitrary keyword arguments.
                Passed to ``mh_next_config``.

        Returns:
            float: The nondimensional end-to-end length.

        Example:
            Calculate the nondimensional end-to-end length under a
            nondimensional force of 23 and ensure that it is within
            one percent of the exact result:

                >>> import numpy as np
                >>> from ufjc import uFJC
                >>> from ufjc.monte_carlo import MHMCMC
                >>> gamma_MHMCMC = MHMCMC(uFJC()). \
                ...     gamma_isotensional_MHMCMC(23, num_processes=1)
                >>> gamma_exact = uFJC().gamma(23, approach='exact')
                >>> np.abs(gamma_exact - gamma_MHMCMC)/gamma_exact < 1e-2
                array([ True])

        """
        gamma_vector = 0
        num_samples = int(kwargs.get('num_samples', int(1e4)))
        for counter in range(1, 1 + num_samples):
            gamma_vector_next = (
                config[-1, :] - config[0, :])/self.N_b
            gamma_vector += (gamma_vector_next - gamma_vector)/counter
            config = self.mh_next_config(config, **kwargs)

        # Return gamma as a scalar
        return la.norm(gamma_vector)

    def gamma_isotensional_MHMCMC(self, eta, **kwargs):
        r"""Parallel Monte Carlo calculation function for
        the isotensional :math:`\gamma(\eta)`.

        Given a nondimensional force, this function computes the
        nondimensional end-to-end length using the average of several
        Monte Carlo calculations run in parallel.
        Given multiple nondimensional forces, this function does the above
        in serial for each nondimensional force.

        Args:
            eta (array_like):
                The nondimensional force(s).
            **kwargs: Arbitrary keyword arguments.
                Passed to ``parallel_calculation``.

        Returns:
            numpy.ndarray: The nondimensional end-to-end length(s).

        Example:
            Calculate the nondimensional end-to-end length under a
            nondimensional force of 23 and ensure that it is within
            one percent of the exact result:

                >>> from ufjc import uFJC
                >>> from ufjc.monte_carlo import MHMCMC
                >>> gamma_MHMCMC = MHMCMC(uFJC()).gamma_isotensional_MHMCMC(23)
                >>> gamma_exact = uFJC().gamma(23, approach='exact')
                >>> np.abs(gamma_exact - gamma_MHMCMC)/gamma_exact < 1e-2
                array([ True])

        """
        eta = self.np_array(eta)
        gamma = np.zeros(eta.shape)
        for index in range(len(eta)):

            # Assume eta is aligned with the z-axis, without loss of generality
            eta_vector = np.array([0, 0, eta[index]])

            # Use the total potential energy Pi (isotensional ensemble)
            self.__bEc = lambda config: self.beta_Pi_config(config, eta_vector)

            # Calculate gamma using an ensemble average
            gamma[index] = self.parallel_calculation(
                self.gamma_isotensional_MHMCMC_serial, **kwargs)

        # Return the resulting averaged gamma
        return gamma
