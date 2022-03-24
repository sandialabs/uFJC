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

        # Collect the files
        tests_dir = path.dirname(__file__)
        files = glob(path.join(tests_dir, './*.py'))
        files += glob(path.join(tests_dir, '../*.py'))
        files += glob(path.join(tests_dir, '../../*.py'))
        files += glob(path.join(tests_dir, '../../docs/*.py'))
        files += glob(path.join(tests_dir, '../examples/*.py'))

        # Check the style in the files
        style = StyleGuide(quiet=False)
        failures = 0
        for file in files:
            print('Testing ' + path.relpath(file) + '........', end='')
            if style.input_file(file):
                failures += 1
            else:
                print('passed!')

        # Check the total number of failues
        if failures > 0:
            print('Failures detected!')
            self.assertEqual(failures, 0)
        else:
            print('All passed!')


if __name__ == '__main__':  # pragma: no cover
    """For command line execution.

    """
    unittest.main()
