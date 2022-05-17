"""An example module computing error as a function of stiffness.

This module contains a ``main`` function that allows the relative L2 error
to be computed for each of the asymptotic approaches, where the
approach to compare against is the exact approach, when available,
and the numerical quadrature approach otherwise.
The results are plotted using ``matplotlib``.
When executed at the command line, flags are passed as keyword arguments.

Example:
    Compute the relative L2 error when using the reduced and
    asymptotic approaches of approximating the isotensional
    single-chain mechanical response of the Morse-FJC model
    as a function of the nondimensional stiffness:

    ::

        python -m ufjc.examples.error --potential 'morse'

"""

from scipy.integrate import quad


def main(**kwargs):
    r"""Main function for the module.

    This is the main function, called when executing the module from the
    command line, also available when importing the module.

    Args:
        **kwargs: Arbitrary keyword arguments.
            Passed to ``uFJC`` instantiation.

    Example:

        .. plot::

            Plot the relative L2 error for a few different potentials:

                >>> import numpy as np
                >>> from ufjc.examples import error
                >>> error.main(potential='harmonic')
                >>> error.main(potential='log-squared')
                >>> error.main(potential='morse')
                >>> error.main(potential='lennard-jones')

    Example:
        Export .csv files for external use:

            >>> from ufjc.examples import error
            >>> error.main(potential='harmonic', csv=1)

    """

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt

    # Import internal modules
    from ufjc import uFJC

    # Enumerate nondimensional energies
    temp = uFJC(**kwargs)
    scale = temp.varepsilon/temp.kappa
    varepsilons = np.logspace(
        np.log(scale*1e1)/np.log(10),
        np.log(scale*1e4)/np.log(10),
        int(kwargs.get('num_points', 3))
    )

    # Allocate space
    kappas = np.zeros(len(varepsilons))
    L_2_rel_error_norm_asymptotic = np.zeros(len(varepsilons))
    L_2_rel_error_norm_reduced = np.zeros(len(varepsilons))

    # Loop over nondimensional energies
    for i, varepsilon in enumerate(varepsilons):

        # Create the single-chain model
        model = uFJC(varepsilon=varepsilon, **kwargs)
        kappas[i] = model.kappa

        # Upper limit of integration
        if hasattr(model, 'lambda_max'):
            upper_lim = model.lambda_max
        else:
            upper_lim = 0.3*model.varepsilon

        # Define the function to compare with
        def gamma_0(eta):
            if model.potential == 'harmonic':
                return model.gamma(eta, approach='exact')
            else:
                return model.gamma(eta, approach='quadrature')

        # Compute the relative L_2 error norm
        L_2_rel_error_norm_asymptotic[i] = np.sqrt(
            quad(
                lambda eta: (
                    model.gamma(eta, approach='asymptotic') - gamma_0(eta)
                )**2, model.minimum_float, upper_lim
            )[0]/quad(
                lambda eta: gamma_0(eta)**2, model.minimum_float, upper_lim
            )[0]
        )
        L_2_rel_error_norm_reduced[i] = np.sqrt(
            quad(
                lambda eta: (
                    model.gamma(eta, approach='reduced') - gamma_0(eta)
                )**2, model.minimum_float, upper_lim
            )[0]/quad(
                lambda eta: gamma_0(eta)**2, model.minimum_float, upper_lim
            )[0]
        )

    # Check if opted to export .csv files
    if (kwargs.get('csv', 0) == 1) is True:
        np.savetxt(
            './fig-error-' + str(model.potential) + '.csv',
            np.vstack((
                L_2_rel_error_norm_asymptotic,
                L_2_rel_error_norm_reduced,
                kappas
            )).T,
            delimiter=","
        )

    # Plot the L_2 error norm if not opted to export .csv files
    else:
        plt.figure()
        plt.loglog(kappas, L_2_rel_error_norm_asymptotic,
                   'o-', label='asymptotic', linewidth=1.5)
        plt.loglog(kappas, L_2_rel_error_norm_reduced,
                   'o-', label='reduced', linewidth=1.5)
        plt.legend()
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$||\gamma - \gamma_0||_2/||\gamma_0||_2$')
        plt.title(model.potential + '-FJC isotensional relative error')
        plt.show()


if __name__ == '__main__':  # pragma: no cover

    # Import external modules
    import argparse

    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--potential', type=str, default='harmonic')

    # Parse arbitrary float args
    for arg in parser.parse_known_args()[1]:
        if arg.startswith('--'):
            parser.add_argument(arg.split('=')[0], type=float)

    # Execute main function with passed args as kwargs
    main(**vars(parser.parse_args()))
