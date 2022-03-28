"""Module-level functionality for running examples.

Example:
    Run the examples in the package:

    ::

        python -m ufjc.examples

    Alternatively:

        >>> from ufjc.examples import __main__

"""

from ufjc.examples import isotensional_approaches
from ufjc.examples import isotensional_asymptotics
isotensional_approaches.main()
isotensional_asymptotics.main(varepsilon_list=[10, 25, 100, 1000])
