# DataAssim

[![Build Status](https://github.com/Alexander-Barth/DataAssim.jl/workflows/CI/badge.svg)](https://github.com/Alexander-Barth/DataAssim.jl/actions)
[![Build Status Windows](https://ci.appveyor.com/api/projects/status/github/Alexander-Barth/DataAssim.jl?branch=master&svg=true)](https://ci.appveyor.com/project/Alexander-Barth/dataassim-jl)
[![codecov](https://codecov.io/github/Alexander-Barth/DataAssim.jl/graph/badge.svg?token=Cwbcb2tnG4)](https://codecov.io/github/Alexander-Barth/DataAssim.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alexander-barth.github.io/DataAssim.jl/stable/)
[![documentation dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexander-barth.github.io/DataAssim.jl/dev/)


The packages implements data assimilation methods:

* (Extended) Kalman Filter
* Incremental 4D-Var
The Julia code is ported from the Matlab/Octave code generated in the frame of the Sangoma project (http://data-assimilation.net/).



## Example

An example of using to package is available as a [jupyter-notebook](https://nbviewer.jupyter.org/github/Alexander-Barth/DataAssim.jl/blob/master/examples/example.ipynb).


## Reference

Most of the algorithms are described in the review article:

Sanita Vetra-Carvalho, Peter Jan van Leeuwen, Lars Nerger, Alexander Barth, M. Umer Altaf, Pierre Brasseur, Paul Kirchgessner, and Jean-Marie Beckers. State-of-the-art stochastic data assimilation methods for high-dimensional non-Gaussian problems. Tellus A: Dynamic Meteorology and Oceanography, 70(1):1445364, 2018. doi: [10.1080/16000870.2018.1445364](https://doi.org/10.1080/16000870.2018.1445364).
