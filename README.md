# DataAssim

[![Build Status](https://github.com/Alexander-Barth/DataAssim.jl/workflows/CI/badge.svg)](https://github.com/Alexander-Barth/DataAssim.jl/actions)
[![Build Status Windows](https://ci.appveyor.com/api/projects/status/github/Alexander-Barth/DataAssim.jl?branch=master&svg=true)](https://ci.appveyor.com/project/Alexander-Barth/dataassim-jl)
[![codecov.io](http://codecov.io/github/Alexander-Barth/DataAssim.jl/coverage.svg?branch=master)](http://codecov.io/github/Alexander-Barth/DataAssim.jl?branch=master)
<!--[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alexander-barth.github.io/DataAssim.jl/stable/)-->
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://alexander-barth.github.io/DataAssim.jl/latest/)


The packages implements various data assimilation methods:

* (Extended) Kalman Filter
* Incremental 4D-Var
* Ensemble Square Root Filter (EnSRF)
* Ensemble Square Root Filter with serial processing of the observations (serialEnSRF)
* Ensemble Transform Kalman Filter (ETKF)
* Ensemble Transform Kalman Filter (EAKF)
* Singular Evolutive Interpolated Kalman ﬁlter (SEIK)
* Error-subspace Transform Kalman Filter (ESTKF)
* Ensemble Kalman Filter (EnKF)

The Julia code is ported from the Matlab/Octave code generated in the frame of the Sangoma project (http://data-assimilation.net/).



## Example

An example of using to package is available as a [jupyter-notebook](https://nbviewer.jupyter.org/github/Alexander-Barth/DataAssim.jl/blob/master/examples/example.ipynb).

https://user-images.githubusercontent.com/9881475/158639340-fe7ae9d7-20ee-4c73-9f02-ade9b89685fa.mp4


## Reference

Most of the algorithms are described in the review article:

Sanita Vetra-Carvalho, Peter Jan van Leeuwen, Lars Nerger, Alexander Barth, M. Umer Altaf, Pierre Brasseur, Paul Kirchgessner, and Jean-Marie Beckers. State-of-the-art stochastic data assimilation methods for high-dimensional non-Gaussian problems. Tellus A: Dynamic Meteorology and Oceanography, 70(1):1445364, 2018. doi: [10.1080/16000870.2018.1445364](https://doi.org/10.1080/16000870.2018.1445364).
