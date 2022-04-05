"""A module for testing  examples within docstrings.

This module checks all other modules in the package
and ensures that the examples within each docstring
are executed properly and obtain the expected output.

"""

# Import external modules
import doctest
import unittest
from os import path
from glob import glob


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

        # Collect the files
        tests_dir = path.dirname(__file__)
        files = glob(path.join(tests_dir, './*.py'))
        files += glob(path.join(tests_dir, '../*.py'))
        if examples:
            files += glob(path.join(tests_dir, '../examples/*.py'))
        else:
            files.remove(glob(path.join(tests_dir, '../swfjc.py'))[0])

        # Check the docstring examples in the files
        failures = 0
        for file in files:
            print('Testing ' + path.relpath(file) + '........', end='')
            if doctest.testfile(file, module_relative=False)[0]:
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
