

function _talagrand_diagram!(nbins,buffer,x,y)
    @assert size(x,1) == size(y,1)
    @assert size(x,2) == size(buffer,1)
    @assert size(nbins,1) == size(buffer,1) + 1

    @inbounds for i = 1:size(x,1)
        buffer .= @view x[i,:]
        sort!(buffer)
        ibin = count(<(y[i]), buffer) + 1
        nbins[ibin] += 1
    end

    return nbins
end

_talagrand_diagram(x,y) = _talagrand_diagram!(zeros(Int,size(x,2)+1),x[1,:],x,y)

"""
    freq = talagrand_diagram(x,y)


Compute the frequency of each bins for a Talagrand Diagram where
`x` is the ensemble (array of the size sample × ensemble size) and `y` are
the obervations (sample).

A U-shaped diagram means that the ensemble spread is too small and that
too many observations fall outside the extremes of the ensemble while a ∩-shaped
diagram mean the spread too larage and too many observations falling near the
center of the ensemble.

## Example

```julia
using PyPlot, DataAssim
Nens = 100
Nsample = 10000
# Ensemble
x = randn(Nsample,Nens)
# Observation
y = randn(Nsample)

freq = DataAssim.talagrand_diagram(x,y)
plot(freq)
plot(ones(size(freq)) / size(freq,1))
ylim(0,ylim()[2])
xlabel("bins")
ylabel("frequency")
```
"""
talagrand_diagram(x,y) = _talagrand_diagram(x,y) / size(x,1)
