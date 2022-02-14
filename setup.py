import setuptools
from distutils.core import setup
setup(
    name='ufjc',
    packages=['ufjc'],
    version='1.0.0',
    description='The uFJC single-chain model Python package.',
    author='Michael R. Buche, Scott J. Grutzik',
    author_email='mrbuche@sandia.gov, sjgrutz@sandia.gov',
    license='BSD',
    url='https://sandialabs.github.io/ufjc',
    keywords=['ufjc', 'polymer', 'single', 'chain', 'freely', 'jointed',
              'model', 'statistical', 'mechanics', 'thermodynamics'],
    install_requires=['numpy', 'scipy'],
    extras_require={
      'docs':['anybadge', 'matplotlib', 'sphinx', 'sphinx-rtd-theme', 'sphinxcontrib-bibtex'],
      'plotting':['matplotlib'],
      'testing':['colorama', 'matplotlib', 'pycodestyle', 'pytest', 'pytest-cov'],},
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    project_urls={
      'Anaconda': 'https://anaconda.org/mrbuche/ufjc',
      'Documentation': 'https://ufjc.readthedocs.io/en/latest/?badge=latest',
      'GitHub': 'https://github.com/sandialabs/ufjc',
      'Issues': 'https://github.com/sandialabs/ufjc/issues',
    },
)
