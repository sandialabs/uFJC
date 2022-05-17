"""An example module involving isotensional asymptotics.

This module contains a ``main`` function that allows the asymptotic
approaches to be compared to exact approaches, when available, and
otherwise to the numerical quadrature approach.
The results are plotted using ``matplotlib``.
When executed at the command line, flags are passed as keyword arguments.

Example:
    Compare the asymptotic approaches to the exact approach for the EFJC
    model for a few nondimensional energy scales:

    ::

        python -m ufjc.examples.asymptotics --varepsilon_list 10 25 100 1000

"""


def main(**kwargs):
    r"""Main function for the module.

    This is the main function, called when executing the module from the
    command line, also available when importing the module.

    Args:
        **kwargs: Arbitrary keyword arguments.
            Passed to ``uFJC`` instantiation.

    Example:

        .. plot::

            Compare the asymptotic approaches for many models and parameters:

                >>> from ufjc.examples import asymptotics
                >>> asymptotics.main(potential='harmonic',
                ...     varepsilon_list=[10, 25, 100, 1000])
                >>> asymptotics.main(potential='log-squared',
                ...     varepsilon_list=[10, 25, 250])
                >>> asymptotics.main(potential='morse',
                ...     varepsilon_list=[10, 25, 250])
                >>> asymptotics.main(potential='lennard-jones',
                ...     varepsilon_list=[1, 2, 5, 15])
                >>> asymptotics.main(potential='mie', n=10, m=4,
                ...     varepsilon_list=[1, 2, 5, 15])
                >>> asymptotics.main(potential='polynomial',
                ...     coefficients=[1, 2, 3],
                ...     varepsilon_list=[10, 25, 100])

    Example:
        Export .csv files for external use:

            >>> from ufjc.examples import asymptotics
            >>> asymptotics.main(potential='harmonic',
            ...     varepsilon_list=[10, 25, 100, 1000], csv=1)
            >>> asymptotics.main(potential='morse',
            ...     varepsilon_list=[10, 25, 100, 1000], csv=1)

    """

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt

    # Import internal modules
    from ufjc import uFJC

    # Option to show plots or export .csv files
    export_csv_not_plot = kwargs.get('csv', 0) == 1
    if export_csv_not_plot is False:
        plt.figure()

    # Repeat for all values of varepsilon
    for varepsilon in kwargs.get('varepsilon_list'):

        # Function to produce a file name for .csv output
        def filename(model, approach):
            return './fig-' + str(model.potential) + '-' + \
                str(approach) + '-varepsilon-' + str(varepsilon) + '.csv'

        # Create the single-chain model
        model = uFJC(varepsilon=varepsilon, **kwargs)

        # If harmonic potentials, use the exact solution
        if model.potential == 'harmonic':
            eta_scale = model.varepsilon
            eta = np.linspace(1e-2, 0.3*model.varepsilon,
                              int(kwargs.get('num_points_eta', 100)))
            if export_csv_not_plot is True:
                data = np.vstack((model.gamma(eta, approach='exact'),
                                  eta/eta_scale)).T
                np.savetxt(filename(model, 'exact'), data, delimiter=",")
            else:
                plt.plot(model.gamma(eta, approach='exact'), eta/eta_scale,
                         label=r'$\varepsilon=$' + str(varepsilon))

        # Otherwise use a numerical calculation
        else:
            if model.potential == 'polynomial':
                eta_scale = model.varepsilon
            else:
                eta_scale = model.eta_max
            eta = np.linspace(1e-2, 0.99*eta_scale,
                              int(kwargs.get('num_points_eta', 100)))
            if export_csv_not_plot is True:
                data = np.vstack((model.gamma(eta, approach='quadrature'),
                                  eta/eta_scale)).T
                np.savetxt(filename(model, 'quad'), data, delimiter=",")
            else:
                plt.plot(model.gamma(eta, approach='quadrature'),
                         eta/eta_scale,
                         label=r'$\varepsilon=$' + str(varepsilon))

        # The asymptotic approaches
        if export_csv_not_plot is True:
            data = np.vstack((model.gamma(eta, approach='asymptotic'),
                              eta/eta_scale)).T
            np.savetxt(filename(model, 'asymptotic'), data, delimiter=",")
        else:
            plt.plot(model.gamma(eta, approach='asymptotic'),
                     eta/eta_scale, 'k:')

        # The reduced asymptotic results
        if export_csv_not_plot is True:
            data = np.vstack((model.gamma(eta, approach='reduced'),
                              eta/eta_scale)).T
            np.savetxt(filename(model, 'reduced'), data, delimiter=",")
        else:
            plt.plot(model.gamma(eta, approach='reduced'),
                     eta/eta_scale, 'k--')

    # Show the results if plotting
    if export_csv_not_plot is False:

        # Create linestyle (approach) legend
        if model.potential == 'harmonic':
            plt.plot(np.nan, np.nan, 'k', label='exact')
        else:
            plt.plot(np.nan, np.nan, 'k', label='quadrature')
        plt.plot(np.nan, np.nan, 'k:', label='asymptotic')
        plt.plot(np.nan, np.nan, 'k--', label='reduced')

        # Plot labels and legend
        plt.legend()
        plt.xlabel(r'$\gamma$')
        plt.ylabel(r'$\eta/\eta_\mathrm{max}$')
        plt.title(model.potential + '-FJC isotensional ' + r'$\gamma(\eta)$')
        plt.show()


if __name__ == '__main__':  # pragma: no cover

    # Import external modules
    import argparse

    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--potential', type=str, default='harmonic')
    parser.add_argument('-l', '--varepsilon_list', nargs='+', type=float,
                        required=True)

    # Parse arbitrary float args
    for arg in parser.parse_known_args()[1]:
        if arg.startswith('--'):
            parser.add_argument(arg.split('=')[0], type=float)

    # Execute main function with passed args as kwargs
    main(**vars(parser.parse_args()))
