"""Module-level functionality for running examples.

Example:
    Run the examples in the package:

    ::

        python -m ufjc.examples

    Alternatively:

        >>> from ufjc.examples import __main__

"""

from ufjc.examples import error
from ufjc.examples import approaches
from ufjc.examples import asymptotics
error.main()
approaches.main()
asymptotics.main(varepsilon_list=[10, 25, 100, 1000])
