"""Module-level functionality for testing installation.

Example:
    Test that the package was installed properly:

    ::

        python -m ufjc.tests

"""

from .test_doc import TestDocstringExamples  # pragma: no cover
TestDocstringExamples().test_docstring_python_examples(examples=False)  # pragma: no cover
