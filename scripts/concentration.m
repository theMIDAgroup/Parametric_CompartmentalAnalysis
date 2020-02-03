% C = concentration(k,A,Ca,ei,t0,C0,t,x,w) computes the solution to
% C'= A*C + k*Ca*ei with ei=[0;..;1;0;..0] and initial condition C(t0)=C0.
%                             i-mo
%..........................................................................
% analytic solution
% C(t) = k * int_t0^t exp(A(t-u)) * Ca(u) * ei du + exp(A(t-t0)) * C0 
%..........................................................................
% k is a scalar
% A is a nc x nc matrix.
% Ca is a function handle accepting a vector as argument and returning a
% vector of the same size.
% ei is a nc x 1 vector of {0,1}.
% t0 is a scalar.
% C0 is a vector of length nc.
% t is the time vector.
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.
%
% C is a matrix of size nc x length(t)
%
% NB: concentration allow directly the computation of the concentration 
% starting from t0 != 0

function C = concentration(k,A,Ca,ei,t0,C0,t,x,w)

    nc = length(A);
    nt = length(t);
    ind=find(ei==1);
    
    C = zeros(nc,nt);

    f = @(u)( k * Ca(u) * getcols(expm((t(1)-u)*A),ind) );
    C(:,1) = expm((t(1)-t0)*A) * C0(:) + quadglv(f,t0,t(1),x,w);

    for n=2:nt 
        
        f = @(u)( k * Ca(u) * getcols(expm((t(n)-u)*A),ind) );
        C(:,n) = expm((t(n)-t(n-1))*A) * C(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
end