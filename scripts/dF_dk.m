% der_k = dF_dk(A,Ca,ei,t0,t,x,w,alpha) computes the derivative with respect
% to k of F(t) = alpha * C(t):
% - C(t) is the analitic solution of the compartmental ODE's system;
% - alpha = ones(1,nc).
%..........................................................................
% ODE's system:
% dC/dt = A*C + k*Ca*ei C(t0)=C0
% analytic solution:
% C(t) = k * int_t0^t exp(A(t-u)) * Ca(u) * ei du + exp(A(t-t0)) * C0 
%..........................................................................
% F(t) = alpha * C(t) = alpha * [ k * int_t0^t exp(A(t-u)) * Ca(u) * ei du + exp(A(t-t0)) * C0 ]
% dF_dk(t) = alpha * dC_dk(t) = alpha * [ int_t0^t exp(A(t-u)) * Ca(u) * ei du ]
%
% A is a nc x nc matrix.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% ei is a nc x 1 vector of {0,1}.
% t0 is a scalar.
% t is the time vector.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
% alpha = ones(1,nc).
%
% dC_dk is a matrix of sixe nc x length(t)
% der_k is a vector of size 1 x length(t)

function der_k = dF_dk(A,Ca,ei,t0,t,x,w,alpha) 

    nc = length(A);
    nt = length(t);
    ind=find(ei==1);
    
    dC_dk = zeros(nc,nt);

    f = @(u)( Ca(u) * getcols(expm((t(1)-u)*A),ind) );
    dC_dk(:,1) = quadglv(f,t0,t(1),x,w);

    for n=2:nt 
        
        f = @(u)( Ca(u) * getcols(expm((t(n)-u)*A),ind) );
        dC_dk(:,n) = expm((t(n)-t(n-1))*A) * dC_dk(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
    der_k = (alpha * dC_dk);
    
end