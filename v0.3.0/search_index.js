var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "DataAssim.jl",
    "title": "DataAssim.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#DataAssim.jl-1",
    "page": "DataAssim.jl",
    "title": "DataAssim.jl",
    "category": "section",
    "text": "Documentation for DataAssim.jl"
},

{
    "location": "#DataAssim.FreeRun",
    "page": "DataAssim.jl",
    "title": "DataAssim.FreeRun",
    "category": "function",
    "text": "x,Hx = FreeRun(ℳ,xi,Q,H,nmax,no)\n\nPerforms a free-run with the model ℳ and nmax time-steps starting at the initial condition xi. Observations at the time steps given in no are  extracted with the observation operator H.\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.KalmanFilter",
    "page": "DataAssim.jl",
    "title": "DataAssim.KalmanFilter",
    "category": "function",
    "text": "x,P = KalmanFilter(xi,Pi,ℳ,Q,yo,R,H,nmax,no)\n\nKalman Filter with the model ℳ and nmax time-steps starting at the initial condition xi and error covariance Pi. Observations yo (and error covariance R) at the time steps given in no are assimilated with the observation operator H.\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.fourDVar",
    "page": "DataAssim.jl",
    "title": "DataAssim.fourDVar",
    "category": "function",
    "text": "x,J = fourDVar(\n        xi,Pi,ℳ,yo,R,H,nmax,no;\n        innerloops = 10,\n        outerloops = 2,\n        tol = 1e-5)\n\nIncremental 4D-Var with the model ℳ and nmax time-steps starting at the initial condition xi and error covariance Pi. Observations yo (and error covariance R) at the time steps given in no are assimilated with the observation operator H.\n\n\n\n\n\n"
},

{
    "location": "#Simulation-driver-1",
    "page": "DataAssim.jl",
    "title": "Simulation driver",
    "category": "section",
    "text": "FreeRun\nKalmanFilter\nfourDVar"
},

{
    "location": "#DataAssim.ETKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.ETKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.ETKF(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.ETKF ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.EnKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.EnKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.EnKF(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.EnKF ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.EnSRF",
    "page": "DataAssim.jl",
    "title": "DataAssim.EnSRF",
    "category": "function",
    "text": "Xa,xa = DataAssim.EnSRF(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.EnSRF ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.EAKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.EAKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.EAKF(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.EAKF ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.SEIK",
    "page": "DataAssim.jl",
    "title": "DataAssim.SEIK",
    "category": "function",
    "text": "Xa,xa = DataAssim.SEIK(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.SEIK ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.ESTKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.ESTKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.ESTKF(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.ESTKF ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.serialEnSRF",
    "page": "DataAssim.jl",
    "title": "DataAssim.serialEnSRF",
    "category": "function",
    "text": "Xa,xa = DataAssim.serialEnSRF(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf and observations y using the DataAssim.serialEnSRF ensemble scheme.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.local_ETKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_ETKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.local_ETKF(Xf,H,y,diagR,part,selectObs,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y using the local DataAssim.ETKF.\n\nInputs:\n\nXf: forecast ensemble (n x N)\nH: observation operator (m x n)\ny: observation (m x 1)\ndiagR: diagonal of the observation error covariance R (m x 1)\npart: vector of integer \"labels\". Every element of the state vector with the same number belong to the same subdomain\nselectObs: callback routine to select observations with a within a subdomain. As input is takes an integer representing the index of the state vector and returns a vector of weights (m x 1). For example:\n\n     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );\n\nor\n\n     selectObs(i) = compact_locfun(L,...\n         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere x and y is the horizontal model grid, xobs and yobs are the locations of the observations and L is a correlation length-scale\n\nOptional inputs:\n\ndisplay: if true, then display progress (false is the default)\nminweight: analysis is performed using observations for which  weights is larger than minweight. (default 1e-8)\nHXf: if non empty, then it is the product H*Xf. In this case, H is not  used\n\nOutput:\n\nXa: the analysis ensemble (n x N)\nxa`: the analysis ensemble mean (n x 1)\n\nSee also: compact_locfun\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.local_EnKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_EnKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.local_EnKF(Xf,H,y,diagR,part,selectObs,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y using the local DataAssim.EnKF.\n\nInputs:\n\nXf: forecast ensemble (n x N)\nH: observation operator (m x n)\ny: observation (m x 1)\ndiagR: diagonal of the observation error covariance R (m x 1)\npart: vector of integer \"labels\". Every element of the state vector with the same number belong to the same subdomain\nselectObs: callback routine to select observations with a within a subdomain. As input is takes an integer representing the index of the state vector and returns a vector of weights (m x 1). For example:\n\n     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );\n\nor\n\n     selectObs(i) = compact_locfun(L,...\n         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere x and y is the horizontal model grid, xobs and yobs are the locations of the observations and L is a correlation length-scale\n\nOptional inputs:\n\ndisplay: if true, then display progress (false is the default)\nminweight: analysis is performed using observations for which  weights is larger than minweight. (default 1e-8)\nHXf: if non empty, then it is the product H*Xf. In this case, H is not  used\n\nOutput:\n\nXa: the analysis ensemble (n x N)\nxa`: the analysis ensemble mean (n x 1)\n\nSee also: compact_locfun\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.local_EnSRF",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_EnSRF",
    "category": "function",
    "text": "Xa,xa = DataAssim.local_EnSRF(Xf,H,y,diagR,part,selectObs,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y using the local DataAssim.EnSRF.\n\nInputs:\n\nXf: forecast ensemble (n x N)\nH: observation operator (m x n)\ny: observation (m x 1)\ndiagR: diagonal of the observation error covariance R (m x 1)\npart: vector of integer \"labels\". Every element of the state vector with the same number belong to the same subdomain\nselectObs: callback routine to select observations with a within a subdomain. As input is takes an integer representing the index of the state vector and returns a vector of weights (m x 1). For example:\n\n     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );\n\nor\n\n     selectObs(i) = compact_locfun(L,...\n         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere x and y is the horizontal model grid, xobs and yobs are the locations of the observations and L is a correlation length-scale\n\nOptional inputs:\n\ndisplay: if true, then display progress (false is the default)\nminweight: analysis is performed using observations for which  weights is larger than minweight. (default 1e-8)\nHXf: if non empty, then it is the product H*Xf. In this case, H is not  used\n\nOutput:\n\nXa: the analysis ensemble (n x N)\nxa`: the analysis ensemble mean (n x 1)\n\nSee also: compact_locfun\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.local_EAKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_EAKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.local_EAKF(Xf,H,y,diagR,part,selectObs,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y using the local DataAssim.EAKF.\n\nInputs:\n\nXf: forecast ensemble (n x N)\nH: observation operator (m x n)\ny: observation (m x 1)\ndiagR: diagonal of the observation error covariance R (m x 1)\npart: vector of integer \"labels\". Every element of the state vector with the same number belong to the same subdomain\nselectObs: callback routine to select observations with a within a subdomain. As input is takes an integer representing the index of the state vector and returns a vector of weights (m x 1). For example:\n\n     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );\n\nor\n\n     selectObs(i) = compact_locfun(L,...\n         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere x and y is the horizontal model grid, xobs and yobs are the locations of the observations and L is a correlation length-scale\n\nOptional inputs:\n\ndisplay: if true, then display progress (false is the default)\nminweight: analysis is performed using observations for which  weights is larger than minweight. (default 1e-8)\nHXf: if non empty, then it is the product H*Xf. In this case, H is not  used\n\nOutput:\n\nXa: the analysis ensemble (n x N)\nxa`: the analysis ensemble mean (n x 1)\n\nSee also: compact_locfun\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.local_SEIK",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_SEIK",
    "category": "function",
    "text": "Xa,xa = DataAssim.local_SEIK(Xf,H,y,diagR,part,selectObs,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y using the local DataAssim.SEIK.\n\nInputs:\n\nXf: forecast ensemble (n x N)\nH: observation operator (m x n)\ny: observation (m x 1)\ndiagR: diagonal of the observation error covariance R (m x 1)\npart: vector of integer \"labels\". Every element of the state vector with the same number belong to the same subdomain\nselectObs: callback routine to select observations with a within a subdomain. As input is takes an integer representing the index of the state vector and returns a vector of weights (m x 1). For example:\n\n     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );\n\nor\n\n     selectObs(i) = compact_locfun(L,...\n         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere x and y is the horizontal model grid, xobs and yobs are the locations of the observations and L is a correlation length-scale\n\nOptional inputs:\n\ndisplay: if true, then display progress (false is the default)\nminweight: analysis is performed using observations for which  weights is larger than minweight. (default 1e-8)\nHXf: if non empty, then it is the product H*Xf. In this case, H is not  used\n\nOutput:\n\nXa: the analysis ensemble (n x N)\nxa`: the analysis ensemble mean (n x 1)\n\nSee also: compact_locfun\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.local_ESTKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_ESTKF",
    "category": "function",
    "text": "Xa,xa = DataAssim.local_ESTKF(Xf,H,y,diagR,part,selectObs,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y using the local DataAssim.ESTKF.\n\nInputs:\n\nXf: forecast ensemble (n x N)\nH: observation operator (m x n)\ny: observation (m x 1)\ndiagR: diagonal of the observation error covariance R (m x 1)\npart: vector of integer \"labels\". Every element of the state vector with the same number belong to the same subdomain\nselectObs: callback routine to select observations with a within a subdomain. As input is takes an integer representing the index of the state vector and returns a vector of weights (m x 1). For example:\n\n     selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );\n\nor\n\n     selectObs(i) = compact_locfun(L,...\n         sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere x and y is the horizontal model grid, xobs and yobs are the locations of the observations and L is a correlation length-scale\n\nOptional inputs:\n\ndisplay: if true, then display progress (false is the default)\nminweight: analysis is performed using observations for which  weights is larger than minweight. (default 1e-8)\nHXf: if non empty, then it is the product H*Xf. In this case, H is not  used\n\nOutput:\n\nXa: the analysis ensemble (n x N)\nxa`: the analysis ensemble mean (n x 1)\n\nSee also: compact_locfun\n\n\n\n\n\n"
},

{
    "location": "#Ensemble-methods-1",
    "page": "DataAssim.jl",
    "title": "Ensemble methods",
    "category": "section",
    "text": "ETKF\nEnKF\nEnSRF\nEAKF\nSEIK\nESTKF\nserialEnSRF\nlocal_ETKF\nlocal_EnKF\nlocal_EnSRF\nlocal_EAKF\nlocal_SEIK\nlocal_ESTKF"
},

{
    "location": "#DataAssim.AbstractModel",
    "page": "DataAssim.jl",
    "title": "DataAssim.AbstractModel",
    "category": "type",
    "text": "Abstract base-class of models. A model should implement forecast step, tangent-linear and adjoint step\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.ModelMatrix",
    "page": "DataAssim.jl",
    "title": "DataAssim.ModelMatrix",
    "category": "type",
    "text": "ℳ = ModelMatrix(M)\n\nLinear model defined by the matrix M.\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.ModelFun",
    "page": "DataAssim.jl",
    "title": "DataAssim.ModelFun",
    "category": "type",
    "text": "ℳ = ModelFun(nonlinear_forecast,tangent_linear_model,adjoint_model)\n\nModel defined by the functions nonlinear_forecast,tangent_linear_model and  adjoint_model.\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.LinShallowWater1DModel",
    "page": "DataAssim.jl",
    "title": "DataAssim.LinShallowWater1DModel",
    "category": "type",
    "text": "ℳ = LinShallowWater1DModel(dt,g,h,L,imax)\n\nLinear 1D shallow water model.\n\nExample\n\ndt = 1.\ng = 9.81\nh = 100\nimax = 101\nL = 10000\nLinShallowWater1DModel(dt,g,h,L,imax)\n\n\n\n\n\n"
},

{
    "location": "#DataAssim.Lorenz63Model",
    "page": "DataAssim.jl",
    "title": "DataAssim.Lorenz63Model",
    "category": "type",
    "text": "ℳ = Lorenz63Model(dt,σ=10.,β = 8/3.,ρ = 28.)\n\nLorenz, 1963 model[1] integrated with a 2nd order Runge-Kutta scheme.\n\n[1] https://doi.org/10.1175/1520-0469(1963)020%3C0130:DNF%3E2.0.CO;2\n\n\n\n\n\n"
},

{
    "location": "#Models-1",
    "page": "DataAssim.jl",
    "title": "Models",
    "category": "section",
    "text": "AbstractModel\nModelMatrix\nModelFun\nLinShallowWater1DModel\nLorenz63Model"
},

{
    "location": "#DataAssim.compact_locfun",
    "page": "DataAssim.jl",
    "title": "DataAssim.compact_locfun",
    "category": "function",
    "text": " fun = compact_locfun(r)\n\nSmooth compact localization function at the (scaled) distance r. fun is zero if r > 2 and one if r is 0. (Gaspari et al. (1999), equation 4.10, [1])\n\n[1] http://dx.doi.org/10.1002/qj.49712555417\n\n\n\n\n\n"
},

{
    "location": "#Utility-functions-1",
    "page": "DataAssim.jl",
    "title": "Utility functions",
    "category": "section",
    "text": "compact_locfun"
},

]}
