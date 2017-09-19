# DataAssim

[![Build Status Linux and macOS](https://travis-ci.org/Alexander-Barth/DataAssim.jl.svg?branch=master)](https://travis-ci.org/Alexander-Barth/DataAssim.jl)
[![Build Status Windows](https://ci.appveyor.com/api/projects/status/github/Alexander-Barth/DataAssim.jl?branch=master&svg=true)](https://ci.appveyor.com/project/Alexander-Barth/ncdatasets-jl)

[![Coverage Status](https://coveralls.io/repos/Alexander-Barth/DataAssim.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Alexander-Barth/DataAssim.jl?branch=master)
[![codecov.io](http://codecov.io/github/Alexander-Barth/DataAssim.jl/coverage.svg?branch=master)](http://codecov.io/github/Alexander-Barth/DataAssim.jl?branch=master)

[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alexander-barth.github.io/DataAssim.jl/stable/)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://alexander-barth.github.io/DataAssim.jl/latest/)

[![DataAssim](http://pkg.julialang.org/badges/DataAssim_0.6.svg)](http://pkg.julialang.org/?pkg=DataAssim)



The packages implements various ensemble Kalman Filter data assimilation methods:

* Ensemble Sqare Root Filter (EnSRF)
* Ensemble Sqare Root Filter with serial processsing of the observations (serialEnSRF)
* Ensemble Transform Kalman Filter (ETKF)
* Ensemble Transform Kalman Filter (EAKF)
* Singular Evolutive Interpolated Kalman ﬁlter (SEIK)
* Error-subspace Transform Kalman Filter (ESTKF)
* Ensemble Kalman Filter (EnKF)

The Julia code is ported from the Matlab/Octave code generated in the frame of the Sangoma project (http://data-assimilation.net/).
