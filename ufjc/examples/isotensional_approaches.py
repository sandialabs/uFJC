"""An example module comparing isotensional approaches.

This module contains a ``main`` function that allows the comparison of
different approaches of obtaining the single-chain mechanical response of
a given single-chain model in the isotensional ensemble.
The results are plotted using ``matplotlib``.
When executed at the command line, flags are passed as keyword arguments.

Example:
    Compare approaches, using 1000 samples for the Monte Carlo approach:

    ::

        python isotensional_approaches --num_samples 1000

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

            Plot and compare each approach for each available potential,
            including examples of polynomial and fully custom potentials,
            for :math:`\varepsilon=88`, keeping the harmonic result (EFJC)
            as a comparison for non-harmonic cases:

                >>> import numpy as np
                >>> from ufjc.examples import isotensional_approaches
                >>> isotensional_approaches.main()
                >>> isotensional_approaches.main(potential='log-squared')
                >>> isotensional_approaches.main(potential='morse')
                >>> isotensional_approaches.main(potential='lennard-jones')
                >>> isotensional_approaches.main(potential='mie', n=10, m=4)
                >>> isotensional_approaches.main(potential='polynomial',
                ...     coefficients=[1, 2, 3])
                >>> isotensional_approaches.main(potential='custom',
                ...     varepsilon=88,
                ...     phi=lambda lambda_: 1 - np.cos(lambda_ - 1),
                ...     eta_link=lambda lambda_: 88*np.sin(lambda_ - 1),
                ...     delta_lambda=lambda eta: np.arcsin(eta/88),
                ...     kappa=88,
                ...     c=2/3,
                ...     lambda_max=1 + np.pi/2,
                ...     eta_max=88)

    """

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt

    # Import internal modules
    from ufjc import uFJC

    # Create the single-chain model
    model = uFJC(**kwargs)

    # Nondimensional forces to apply
    if kwargs.get('eta_max', None) is not None:
        eta_max = kwargs.get('eta_max', None)
    elif hasattr(model, 'eta_max'):
        eta_max = 0.95*model.eta_max
    else:
        eta_max = (0.5 + model.kappa/4)
    eta = np.linspace(1e-2, eta_max, 100)
    eta_mc = np.linspace(eta[0], eta[-1], int(kwargs.get('num_forces', 5)))

    # Plot each available approach for the isotensional gamma(eta)
    plt.figure()
    if model.potential == 'harmonic':
        plt.plot(model.gamma(eta, approach='exact'), eta,
                 label='exact', linewidth=2.5)
    else:
        model_h = uFJC(varepsilon=model.kappa, potential='harmonic')
        gamma_h = model_h.gamma(eta, approach='exact')
        plt.plot(gamma_h, eta, 'k-', label='(EFJC exact)', linewidth=1.5)
    plt.plot(model.gamma(eta, approach='asymptotic'), eta,
             '--', label='asymptotic', linewidth=2.5)
    plt.plot(model.gamma(eta, approach='reduced'), eta,
             '--', label='reduced', linewidth=2.5)
    plt.plot(model.gamma(eta, approach='quadrature'), eta,
             ':', label='quadrature', linewidth=4)
    plt.plot(model.gamma(eta_mc, approach='monte carlo', **kwargs), eta_mc,
             'ko', label='Monte Carlo')
    plt.legend()
    plt.xlabel(r'$\gamma$')
    plt.ylabel(r'$\eta$')
    plt.xlim([-0.03, 1.08*model.gamma(eta[-1], approach='asymptotic')])
    plt.title(model.potential + '-FJC isotensional ' + r'$\gamma(\eta)$')
    plt.show()


if __name__ == '__main__':
    """For command line execution.

    """
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
