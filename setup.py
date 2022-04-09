from distutils.core import setup
from setuptools import find_packages


def read(fname):
    with open(fname) as fp:
        content = fp.read()
    return content


setup(
    name='ufjc',
    version='1.3.0',
    package_dir={'ufjc': 'ufjc'},
    packages=find_packages(),
    description='The Python package for the uFJC single-chain model.',
    long_description=read("README.rst"),
    author='Michael R. Buche, Scott J. Grutzik',
    author_email='mrbuche@sandia.gov, sjgrutz@sandia.gov',
    url='https://sandialabs.github.io/ufjc',
    license='BSD-3-Clause',
    keywords=['ufjc', 'polymers', 'statistical mechanics', 'thermodynamics'],
    install_requires=['numpy', 'pycodestyle', 'scipy'],
    extras_require={
      'docs': ['matplotlib',
               'sphinx', 'sphinx-rtd-theme', 'sphinxcontrib-bibtex'],
      'plotting': ['matplotlib'],
      'testing': ['matplotlib', 'pytest', 'pytest-cov'],
      'all': ['matplotlib', 'pytest', 'pytest-cov',
              'sphinx', 'sphinx-rtd-theme', 'sphinxcontrib-bibtex']
    },
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    project_urls={
      'Anaconda': 'https://anaconda.org/mrbuche/ufjc',
      'Documentation': 'https://ufjc.readthedocs.io',
      'GitHub': 'https://github.com/sandialabs/ufjc',
      'Issues': 'https://github.com/sandialabs/ufjc/issues',
    },
)
