% der_A = dF_dA(k,A,Ca,ei,t0,C0,t,x,w,alpha) computes the derivative 
% with respect to A of F(t) = alpha * C(t):
% - C(t) is the analitic solution of the compartmental ODE's system;
% - alpha = ones(1,nc).
%..........................................................................
% ODE's system:
% dC/dt = A*C + k*Ca*ei C(t0)=C0
% analytic solution:
% C(t) = k * int_t0^t exp(A(t-u)) * Ca(u) * ei du + exp(A(t-t0)) * C0 
%..........................................................................
% F(t) = alpha * C(t) = alpha * [ k * int_t0^t exp(A(t-u)) * Ca(u) * ei du + exp(A(t-t0)) * C0 ]
% dF_dA(t) = alpha * [(dC/dA(t)).H] = alpha * [int_t0^t kron(C(u).',exp((t-u)A) du] * H(:)
%
% k is a scalar.
% A is a nc x nc matrix.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% ei is a nc x 1 vector of {0,1}.
% t0 is a scalar.
% C0 is a vector of length nc.
% t is the time vector.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
% alpha = ones(1,nc).
%
% dC_dA is an array of size nc x nc^2 x length(t).
% Each nc x nc^2 matrix gives the values of the derivative at times t.
% der_A is a matrix of size nc^2 x length(t).

function der_A = dF_dA(k,A,Ca,ei,t0,C0,t,x,w,alpha)
  
    nc = length(A);
    nt = length(t);
   
    dC_dA = zeros(nc,nc^2,nt);
        
    C = concentration(k,A,Ca,ei,t0,C0,t,x,w);
    
    Cu = @(u)( concentration(k,A,Ca,ei,t0,C0,u,x,w) );
    f = @(u)( kron((Cu(u)).',expm((t(1)-u)*A)) );
    
    dC_dA(:,:,1) = quadglv(f,t0,t(1),x,w);
    
    for n=2:nt;
        
        Cu = @(u)( concentration(k,A,Ca,ei,t(n-1),C(:,n-1),u,x,w) );
        f = @(u)( kron((Cu(u)).',expm((t(n)-u)*A)) );
        
        dC_dA(:,:,n) = expm((t(n)-t(n-1))*A) * dC_dA(:,:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
   
    der_A = reshape(alpha*reshape(dC_dA,nc,nc^2*nt),nc^2,nt);
     
end    