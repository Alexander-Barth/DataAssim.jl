using Statistics, Random, Test

# Algorithm 1: Algorithm of probabilistic resampling (PR)

function PR(w,Ñₑ)
    Nₑ = length(w)
    ŵ = similar(w, eltype(w), (Nₑ,))
    ŵ[1] = w[1]

    for j = 2:Nₑ
        ŵ[j] = ŵ[j-1] + w[j]
    end

    k = zeros(Int,Ñₑ)
    c = 1
    for j = 1:Ñₑ
        u = rand()

        while u > ŵ[c]
            c = c+1
        end

        k[j] = c
        c = 1
    end
    return k
end

k = PR([0,0.1,0.9],10000);
@test sum(k .== 1) == 0

# Algorithm 2 Algorithm of stochastic universal resampling (SUR)
function SUR(w,Ñₑ)
    Nₑ = length(w)

    ŵ = similar(w, eltype(w), (Nₑ,))
    ŵ[1] = w[1]

    for j = 2:Nₑ
        ŵ[j] = ŵ[j-1] + w[j]
    end

    k = zeros(Int,Ñₑ)
    # generate a random number
    u = rand()/Nₑ
    c = 1
    @show ŵ

    for j = 1:Ñₑ
        @show c
        while u > ŵ[c]
            @show c,u,ŵ[c]
            c = c+1
        end

        k[j] = c
        u = u + 1/Nₑ
        c = 1
    end

    return k
end

#k = SUR([0,0.1,0.9],10000);
