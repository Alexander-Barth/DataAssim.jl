
"""
    FreeRun!(ℳ, xi, Q, 𝓗, nmax, no, x, Hx)

Performs forward integration of the model `ℳ` starting from initial state `xi` up to time step `nmax`.

- `x[:,n] = ℳ(n, x[:,n])` is the state estimated by the model at time n.
- `𝓗(obsindex, x)` mapping.
- `nmax` is total number of integration of the model 
- `no` is a vector of time indices where observations exist.
- `x` is the matrix to store the state trajectory.
- `Hx` is the vector to store the simulated observations at times `no`.

Modifies `x` and `Hx` in place.
"""

function FreeRun!(ℳ,xi,Q,𝓗::AbstractModel,nmax,no,x,Hx)
    obsindex = 1;

    for n=1:nmax+1
        if n == 1
            x[:,1] = xi;
        else
            x[:,n] = ℳ(n-1,x[:,n-1]);
        end

        if obsindex <= length(no) && n == no[obsindex]
            # extract observations
            Hx[obsindex] = 𝓗(obsindex,x[:,n])
            obsindex = obsindex +1;
        end
    end
end

"""
    x,Hx = FreeRun(ℳ,xi,Q,H,nmax,no)

Performs a free-run with the model `ℳ` and `nmax` time-steps starting at the
initial condition `xi`. Observations at the time steps given in `no` are
extracted with the observation operator `H`.
"""
function FreeRun(ℳ,xi,Q,𝓗,nmax,no)
    T = eltype(xi)
    x = zeros(T,size(xi,1),nmax+1);
    Hx = Vector{Vector{T}}(undef,length(no))
    FreeRun!(ℳ,xi,Q,𝓗,nmax,no,x,Hx)
    return x,Hx
end

"""
    J = costfun(xi, Pi, ℳ, xa, yo, R, 𝓗, nmax, no, x, Hx)

Computes the 4D-Var cost function:

    J = (xa - xi)' * Pi⁻¹ * (xa - xi) + Σ (yo[i] - 𝓗(x[no[i]]))' * R[i]⁻¹ * (yo[i] - 𝓗(x[no[i]]))

Arguments:
- `xi`: background state (initial state)
- `Pi`: background error covariance
- `xa`: current state estimate
- `yo`: list of observations
- `R`: list of observation error covariance matrices
- `𝓗`: observation operator (Maping from state space to observation space)
- `x`, `Hx`: temporary storage for model trajectory and simulated observations

Returns:
- `J`: scalar cost function value
"""
function costfun(xi,Pi,ℳ,xa,yo,R,𝓗,nmax,no,x,Hx)
    FreeRun!(ℳ,xi,[],𝓗,nmax,no,x,Hx);
    # cost function
    tmp = x[:,1] - xa;
    J = tmp' * (Pi \ tmp);
    for i = 1:length(no)
        tmp = yo[i] - Hx[i];
        J = J + tmp' * (R[i] \ tmp);
    end
    return J
end


function gradient(xi,dx0,x,Pi,ℳ,yo,R,𝓗,nmax,no)

    dx = zeros(size(xi,1),nmax+1);
    dx[:,1] = dx0;
    obsindex = length(no);

    for n=1:nmax
        dx[:,n+1] = tgl(ℳ,n,x[:,n],dx[:,n]);
    end

    lambda = zeros(size(xi,1),nmax+2);
    for n=nmax+1:-1:1
        lambda[:,n] = adj(ℳ,n,x[:,n],lambda[:,n+1]);

        if obsindex > 0 && n == no[obsindex]
            lambda[:,n] = lambda[:,n] +
                adj(𝓗,n,x[:,n],R[obsindex] \ (yo[obsindex] - tgl(𝓗,n,x[:,n],dx[:,n]+x[:,n])))
            obsindex = obsindex - 1;
        end
    end

    grad = 2 * (Pi \ ((dx[:,1]+x[:,1]) - xi)) - 2 * lambda[:,1];

    return grad,lambda
end

function Jacobien(xi,dx,x,Pi,ℳ,R,𝓗,nmax,no)
    Jac = []
    obsindex = length(no)
    push!(Jac, cholesky(Pi).U \ dx)
    for n=nmax+1:-1:1
        if obsindex > 0 && n == no[obsindex]
            push!(Jac, cholesky(R[obsindex]).U \ tgl(𝓗,n,x[:,n], tgl(ℳ,n,x[:,n],dx)))
            obsindex = obsindex - 1;
        end
    end
    result = vcat(Jac...)
    return result
end

function Adj_Jacobien(xi,dx,x,Pi,ℳ,R,𝓗,nmax,no)
    Jac = zeros(size(xi,1))
    obsindex = length(no)
    index = size(Pi,1)
    Jac = (cholesky(Pi).U)' \ dx[1:index]
    for n=nmax+1:-1:1
        if obsindex > 0 && n == no[obsindex]
            m = size(R[obsindex],1)
            Jac += adj(ℳ,n,x[:,n], adj(𝓗,n,x[:,n],(cholesky(R[obsindex]).U)' \ dx[index+1:index+m]))
            index +=m 
            obsindex = obsindex - 1;
        end
    end 
    return Jac
end



"""
    compute_eigenvalues(A_mul!, n)

Calcule et affiche la matrice et les valeurs propres de l'opérateur défini
par `A_mul!` de dimension `n`.

Arguments:
- `A_mul!`: fonction mutating telle que A_mul!(y, x) = A*x
- `n`: dimension de l'opérateur
"""
function compute_eigenvalues(A_mul!, n)
    # Initialiser la matrice pleine
    A = zeros(n, n)
    
    # Vecteurs de travail
    x = zeros(n)
    y = zeros(n)
    
    # Construire chaque colonne en appliquant A_mul! aux vecteurs canoniques
    for j in 1:n
        x .= 0.0
        x[j] = 1.0
        A_mul!(y, x)
        A[:, j] .= y
    end
    
    # Affichage de A
    display(A)
    println(A == diag(A))

    # Calcul des valeurs propres
    eigvals = eigen(A).values
    display(eigvals)
    # Graphe des valeurs propres
    #p = Plots.plot(eigvals)
    #Plots.display(p)

end

"""
    x,J = fourDVar(
            xi,Pi,ℳ,yo,R,H,nmax,no;
            innerloops = 10,
            outerloops = 2,
            tol = 1e-5)

Incremental 4D-Var with the model `ℳ` (`AbstractModel`) and `nmax` time-steps starting at the
initial condition `xi` and error covariance `Pi` with the specified numbers of inner
and outer loops.
Observations `yo` (vector of vectors) and error covariance `R` (vector of matrices) at the time steps given in `no` are
assimilated with the observation operator `H` (`AbstractModel`).
"""
function fourDVar(
    xi::AbstractVector,Pi,ℳ,yo::AbstractVector,R::AbstractVector,𝓗,nmax,no;
    outerloops = 1,
    tol = 1e-8)

    xa = float(xi)
    T = eltype(xa)
    x = zeros(size(xi,1),nmax+1);
    Hx = Vector{Vector{T}}(undef,length(no))

    J = zeros(outerloops)
    b = zeros(size(xi))

    for i=1:outerloops

        J[i] = costfun(xi,Pi,ℳ,xa,yo,R,𝓗,nmax,no,x,Hx);

        # dx increment relative to xi

        grad(dx) = gradient(xi,dx,x,Pi,ℳ,yo,R,𝓗,nmax,no)[1];
        b .= grad(zeros(size(xi)));

        Jacob(dx) = Jacobien(xi,dx,x,Pi,ℳ,R,𝓗,nmax,no)
        AdjJacob(dx) = Adj_Jacobien(xi,dx,x,Pi,ℳ,R,𝓗,nmax,no)

        #function A_mul!(y, dx)
        #    y .= grad(dx) - b
        #    y
        #end

        function A_mul!(y, dx)
            JJ = zeros(eltype(Jacob(dx)))
            JJ = Jacob(dx)
            y .= AdjJacob(JJ)
            y
        end
        
        Hess = LinearOperators.LinearOperator(Float64, length(b), length(b), true, false, A_mul!)
        
        
        global Aerrors = [] 

        function A_error(workspace::KrylovWorkspace, Hess,b)
                error = zeros(length(b))
                (xast,_) = Krylov.cg(Hess, -b)
                error = xast - workspace.x
                mul!(error, Hess, error)
                error = error'*(xast - workspace.x)
                push!(Aerrors, error)
            return false
        end

        cg_callback(workspace) = A_error(workspace, Hess,b)



        
        
        (dxa, stats) = Krylov.cg(Hess, -b, 
                atol=tol, rtol=tol, 
                history=true, 
                callback = cg_callback
                )
        # add increment to dxa
        xa .= xa + dxa;

        
        compute_eigenvalues(A_mul!, length(b))
    end

    return xa,J,Aerrors
end
