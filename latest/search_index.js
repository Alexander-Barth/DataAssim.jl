var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "DataAssim.jl",
    "title": "DataAssim.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DataAssim.jl-1",
    "page": "DataAssim.jl",
    "title": "DataAssim.jl",
    "category": "section",
    "text": "Documentation for DataAssim.jl"
},

{
    "location": "index.html#DataAssim.ETKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.ETKF",
    "category": "function",
    "text": "    Xa,xa = \"ensemble_analysis\"(Xf,HXf,y,R,H,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf\nand observations y using various ensemble scheme.\nThe function name \"ensemble_analysis\" can be EnSRF, EAKF, ETKF, SEIK, ESTKF, serialEnSRF or EnKF.\n\nInput arguments:\n\nXf: forecast ensemble (n x N)\nHXf: the observation operator applied on the ensemble (product H*Xf)\ny: observations (m)\nR: observation error covariance  (m x m).\nH: operator (m x n). Except for the serialEnSRF it is never used and can be empty\n\nOptional keywords arguments:\n\ndebug: set to true to enable debugging. Default (false) is no debugging.\ntolerance: expected rounding error (default 1e-10) for debugging checks. This is not used if debug is false.\n\nOutput arguments:\n\nXa: the analysis ensemble (n x N)\nxa: the analysis ensemble mean (n)\n\nNotations follows: Sangoma D3.1 http://data-assimilation.net/Documents/sangomaDL3.1.pdf\n\n\n\n"
},

{
    "location": "index.html#DataAssim.local_ETKF",
    "page": "DataAssim.jl",
    "title": "DataAssim.local_ETKF",
    "category": "function",
    "text": "[Xa,xa] = local_ensemble_analysis(...    Xf,H,y,diagR,part,selectObs,method,...)\n\nComputes analysis ensemble Xa based on forecast ensemble Xf using the observation y.\n\nInputs: Xf: forecast ensemble (n x N) H: observation operator (m x n) y: observation (m x 1) diagR: diagonal of the observation error covariance R (m x 1) part: vector of integer \"labels\". Every element of the state vector with the   same number belong to the same subdomain selectObs: callback routine to select observations with a within a subdomain.   As input is takes an integer representing the index of the state vector and   returns a vector of weights (m x 1).   For example:      selectObs(i) = exp(- ((x[i] - xobs[:]).^2 + (y(i) - yobs[:]).^2)/L^2 );   or      selectObs(i) = sangoma_compact_locfun(L,...          sqrt((x[i] - xobs[:]).^2 + (y[i] - yobs[:]).^2));\n\nwhere:      x and y is the horizontal model grid      xobs and yobs are the localtion of the observations      L is a correlation length-scale\n\nmethod: method is one analysis schemes implemented sangoma_ensemble_analysis   (except for EnSRF)\n\nOptional inputs: \'display\', display: if true, then display progress (false is the default) \'minweight\', minweight: analysis is performed using observations for which    weights is larger than minweight. (default 1e-8) \'HXf\', HXf: if non empty, then it is the product H Xf. In this case, H is not    used\n\nOutput: Xa: the analysis ensemble (n x N) xa: the analysis ensemble mean (n x 1)\n\nSee also: ensemble_analysis, sangoma_compact_locfun\n\n\n\n"
},

{
    "location": "index.html#Datasets-1",
    "page": "DataAssim.jl",
    "title": "Datasets",
    "category": "section",
    "text": "ETKF\nlocal_ETKF"
},

]}
