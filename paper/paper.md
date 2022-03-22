---
title: 'ufjc: The Python package for the uFJC single-chain model'
tags:
  - python
  - statistical mechanics
  - thermodynamics
  - polymers
authors:
  - name: Michael R. Buche^[mrbuche\@sandia.gov]
    orcid: 0000-0003-1892-0502
    affiliation: 1
  - name: Scott J. Grutzik
    orcid: 0000-0002-6490-3941
    affiliation: 1
affiliations:
 - name: Materials and Failure Modeling, Sandia National Laboratories, Albuquerque, NM 87185, USA
   index: 1
date: 16 March 2022
bibliography: paper.bib
---

# Summary

``ufjc`` is the Python package that implements the $u$FJC model for single polymer chains.
The $u$FJC model replaces the rigid links of the freely-jointed chain (FJC) model [@rubinstein2003polymer] with flexible links of potential energy function $u$ [@buche2021chain].
This replacement allows the stretching of polymer chains to include bond stretching, and depending on the potential choice, bond breaking [@buche2022on].
Through a robust implementation the $u$FJC single-chain model, the ``ufjc`` package 
(1) allows the stretching of polymer chains to be efficiently modeled,
(2) is readily integrated into single-chain-based polymer network constitutive models, and
(3) provides the tools to study the fundamental statistical thermodynamics of the model system.
``ufjc`` utilizes efficient capabilities provided by ``numpy`` [@numpy] and ``scipy`` [@scipy]; visualization using ``matplotlib`` [@matplotlib] is recommended.

# Basic usage

The model class is first imported from the package,

```python
>>> from ufjc import uFJC
```

and then a model instance may be created,

```python
>>> model = uFJC(N_b=8, potential='morse', varepsilon=88)
```

Here, ``model`` is an instance of the $u$FJC single-chain model with 8 links, where each link is assigned the @morse1929diatomic potential with a nondimensional energy scale of 88 [@buche2022on].
The model instance has methods corresponding to thermodynamic functions, and results are returned as a ``numpy.ndarray``.
For example, the nondimensional end-to-end length (extension) of the chain, $\gamma$, is computed when a nondimensional force of $\eta=0.55$ is applied at the ends of the chain:

```python
>>> model.gamma(0.55)
array([0.18780929])
```

Other methods include the nondimensional free energy of the chain as a function of extension, the nondimensional probability functions of chain extensions at equilibrium, and the net rate of breaking a chain as a function of extension.
Optional keyword arguments in each method are available in order to specify certain features, such as the calculation approach or the thermodynamic ensemble.
For example, evaluate the nondimensional equilibrium radial distribution function using the reduced asymptotic approach:

```python
>>> model.nondim_g_eq([0, 0.23, 0.88], approach='reduced')
array([0.00000000e+00, 2.64707367e+00, 1.24115933e-04])
```

These functions are plotted over an exemplary range of parameters in \autoref{fig:example}.

\newpage

![
  Thermodynamic functions for the Morse-FJC model with base parameters $N_b=8$, $\alpha=1$, and $\varepsilon=88$.
  (top) The nondimensional single-chain mechanical response (force versus extension) for an increasing nondimensional potential energy scale.
  An x denotes the maximum force that the chain can sustain, where additional force would cause breaking.
  (bottom) The nondimensional equilibrium radial distribution function versus extension for an increasing number of links.
\label{fig:example}](figure.png)

# Statement of need

``ufjc`` was developed for researchers to effectively model single polymer chains and networks of polymer chains.
Researchers have historically relied on the freely-jointed chain (FJC) model, which captures fundamental physics while maintaining analytic simplicity.
The FJC model consists of a number of rigid phantom links connected in series by penalty-free hinges; the increasing force required to extend the chain is due to entropy reduction, and is described by the Langevin function [@rubinstein2003polymer].
Researchers have become interested in modeling polymers up to and including failure, but unfortunately the rigid links of the FJC model are neither capable of modeling bond stretching nor breaking.
The rigid links can be made flexible according to some potential energy function $u$, but this causes the model to become analytically intractable, except for the special case of a harmonic potential and the isotensional ensemble [@balabaev2009extension; @manca2012elasticity; @buche2020statistical].
Quite recently, an asymptotically-correct statistical thermodynamic theory [@buche2021fundamental] has been successfully applied to this $u$FJC model, resulting in accurate, analytically-tractable relations for the model for steep link potentials [@buche2022on].
The ``ufjc`` package robustly implements this approach, enabling the results presented by @buche2022on while providing many additional functionalities.
The object-oriented structure of ``ufjc`` allows it to be readily integrated into other packages as a sub-package, such as polymer network constitutive models that build up from a single-chain models.
Correspondingly, ``ufjc`` will streamline ongoing constitutive model development, and improve existing constitutive models through reimplementation [@buche2021chain].

\newpage

# Acknowledgements

This work was supported by the Laboratory Directed Research and Development (LDRD) program at Sandia National Laboratories under project 222398.
Sandia National Laboratories is a multi-mission laboratory managed and operated by National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of Honeywell International, Inc., for the U.S. Department of Energy's National Nuclear Security Administration under contract DE-NA0003525.
This paper describes objective technical results and analysis.
Any subjective views or opinions that might be expressed in the paper do not necessarily represent the views of the U.S. Department of Energy or the United States Government.

# References
