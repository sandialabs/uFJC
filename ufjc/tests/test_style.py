"""A module for testing the Python style.

This module checks all other modules in the package (including itself)
for adherence to the PEP 8 style standard for Python code.

"""

# Import external modules
import unittest
from os import path
from glob import glob
from pycodestyle import StyleGuide


class TestCodeStyle(unittest.TestCase):
    """Class to test conformance of code to PEP style standards.

    """
    def test_pep8_conformance(self):
        """Function to test conformance of code to PEP style standards.

        """
        tests_dir = path.dirname(__file__)
        files = glob(path.join(tests_dir, './*.py'))
        files += glob(path.join(tests_dir, '../*.py'))
        files += glob(path.join(tests_dir, '../../*.py'))
        files += glob(path.join(tests_dir, '../../docs/*.py'))
        files += glob(path.join(tests_dir, '../examples/*.py'))
        style = StyleGuide(quiet=False)
        for file in files:
            self.assertFalse(style.input_file(file))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
