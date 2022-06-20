struct GaussianSqrt
    mean
    SqrtCovar
end

rand(g::GaussianSqrt) = g.mean + g.SqrtCovar * randn(size(g.SqrtCovar,2))


struct Gaussian
    mean
    covar
end

GaussianSqrt(g::Gaussian) = GaussianSqrt(g.mean,cholesky(g.covar).U)
rand(g::Gaussian) = rand(GaussianSqrt(g))
