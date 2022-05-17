"""A module for testing examples within docstrings.

This module checks all other modules in the package
and ensures that the examples within each docstring
are executed properly and obtain the expected output.

"""

# Import external modules
import doctest
import unittest
from os import path
from glob import glob
from sys import platform


class TestDocstringExamples(unittest.TestCase):
    """Class to test examples within docstrings.

    """
    def test_docstring_python_examples(self, examples=True):
        """Function to test examples within docstrings.

       Args:
           examples (bool, optional, default=True): Whether to test examples.

        Returns:
            int: The total number of failures from all files.

        """
        tests_dir = path.dirname(__file__)
        files = glob(path.join(tests_dir, './*.py'))
        files += glob(path.join(tests_dir, '../*.py'))
        if examples:
            files += glob(path.join(tests_dir, '../examples/*.py'))
        else:  # pragma: no cover
            if platform == 'win32' or platform == 'win64':
                files.remove(glob(path.join(tests_dir, '..\\swfjc.py'))[0])
            else:
                files.remove(glob(path.join(tests_dir, '../swfjc.py'))[0])
        for file in files:
            self.assertFalse(doctest.testfile(file, module_relative=False)[0])


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
