"""Module-level functionality for testing installation.

Example:
    Test that the package was installed properly:

    ::

        python -m ufjc.tests

    Alternatively:

        >>> from ufjc.tests import __main__

"""

from .test_doc import TestDocstringExamples
TestDocstringExamples().test_docstring_python_examples(examples=False)
