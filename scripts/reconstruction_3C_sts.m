function [k1x,k2x,k3x,k4x,k5x,k6x,k7x,relerr,nit,Cxdata] = reconstruction_3C_sts(Cdata,Ca,t)

% 3-COMPARTMENT C1,C2,C3 COMPARTMENTAL SYSTEM with SINGLE INPUT FUNCTION Ca
% KIDNEYS COMPATMENTAL MODEL               
%                              Ca(=arterial IF)
%                               |
%                            K2 |
%                               V
%       K1    ______    K4    ______          ______
%      --->  |      |  --->  |      |   K6   |      |   K7
%  Ca   K3   |  C1  |   K5   |  C2  |  --->  |  C3  |  --->  Cu(=bladder)
%      <---  |______|  <---  |______|        |______|
%             C1=free      C2=metabolized   C3=tubular          

% The ODE system is C'=A*C+[K1*Ca;K2*Ca;0] with initial condition C(0)=[0;0;0]
% The solution is C(t)=K1*int( expm((t-u)A)*Ca(u) , u=t0..t) [e1] + K2*int( expm((t-u)A)*Ca(u) , u=t0..t) [e2]

%% Inizial guessed parameters for Gauss-Newton's algorithm

% Fraction of blood volume
Vb = 0.2;
alpha = [1-Vb,1-Vb,1-Vb]; % --> need to sum different compartment C_.. and to obtain C_tot

% Initial values
t0 = 0; C0 = [0;0;0];

% Initial guess
attempt = 1;
[k1x,k2x,k3x,k4x,k5x,k6x,k7x] = deal(2e-1,1e-2,5e-2,2,4,1e-2,1e-4);
%(2e-1*rand(1),1e-2*rand(1),5e-2*rand(1),1*rand(1),1*rand(1),1e-2*rand(1),1e-4*rand(1)); % K6=10^2*K7
% (2e-1,1e-2,5e-2,2,4,1e-2,1e-4)
x = [k1x;k2x;k3x;k4x;k5x;k6x;k7x];

% Gauss-Legendre quadrature
ngl = 4; [glnodes,glweights] = gauss_legendre(ngl);

% Solve DIRECT PROBLEM (with initial parameters)
Ax = [[-(k3x+k4x);k4x;0],[k5x;-(k5x+k6x);k6x],[0;0;-k7x]];

Cx1 = concentration(k1x,Ax,Ca,[1;0;0],t0,C0,t,glnodes,glweights); %[Ca;0;0] on e1
Cx2 = concentration(k2x,Ax,Ca,[0;1;0],t0,C0,t,glnodes,glweights); %[0;Ca;0] on e2
Cx = Cx1+Cx2;
Cxdata = (alpha*Cx + Vb*Ca(t)).';

% Relative error
relerr = norm(Cdata-Cxdata)/norm(Cdata);

%% Regularized GAUSS - NEWTON method 

% Iteration numbers
nit = 0; 
iter = 0;
nit_max = 30;

% Lower bound for the relative error
toll = 8*1e-2; % sqrt(sum(Cdata))/norm(Cdata);

% Initializations
vect_relerr = zeros(2*nit_max,1);
cell_K = cell(1,2*nit_max);

while relerr>toll || nit==0 % <-- I want the algorithm to make at least one iteration! 
    
    % Compute the matrix of the derivative of dot(alpha,C) with respect to K and A
    % Derivative with respect to k1
    aux1 = dF_dk(Ax,Ca,[1;0;0],t0,t,glnodes,glweights,alpha).';
    % Derivative with respect to k2
    aux2 = dF_dk(Ax,Ca,[0;1;0],t0,t,glnodes,glweights,alpha).';
    % Derivative with respect to A
    aux3 = dF_dA(k1x,Ax,Ca,[1;0;0],t0,C0,t,glnodes,glweights,alpha).';
    aux4 = dF_dA(k2x,Ax,Ca,[0;1;0],t0,C0,t,glnodes,glweights,alpha).';
    M = [aux1 aux2 aux3+aux4];
    M = [M(:,1),M(:,2),-M(:,3),M(:,4)-M(:,3),M(:,6)-M(:,7),M(:,8)-M(:,7),-M(:,11)];
    % M=[dalphaCdK1,dalphaCdK2,dalphaCdK3,dalphaCdK4,dalphaCdK5,dalphaCdK6,dalphaCdK7]
    % considering Ax=[[-(K3x+K4x);K4x;0],[K5x;-(K5x+K6x);K6x],[0;0;-K7x]];
        
    % Regularization parameter --------------------------------------------   
    vect_lambda = 1e4:5e2:1e6; %values for lambda to compute V_lambda
    vect_lambda = vect_lambda';
    
    % Generalized Cross Validation
    [r,~] = GCV(M,Cdata-Cxdata,vect_lambda);
    %----------------------------------------------------------------------
    % '\' solve the linear system Mh=D with D=Cdata-Cxdata h=increment
    % ==> Tychonov regularization: (rI+M'*M)h=M'*D
    h = (r*diag([1,1,1,1,1,1,1])+M.'*M)\(M.'*(Cdata-Cxdata));
    
    % Check for positive coefficients
    xph = x+h;
    if any(xph<=0)
        h(xph<=0)=0;
    end
    x = x+h;
    
    % Refresh the parameters
    k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4);
    k5x = x(5); k6x = x(6); k7x = x(7);
    
    % Solve DIRECT PROBLEM --> refresh data
    Ax = [[-(k3x+k4x);k4x;0],[k5x;-(k5x+k6x);k6x],[0;0;-k7x]];
    Cx1 = concentration(k1x,Ax,Ca,[1;0;0],t0,C0,t,glnodes,glweights); %[Ca;0;0] on e1
    Cx2 = concentration(k2x,Ax,Ca,[0;1;0],t0,C0,t,glnodes,glweights); %[0;Ca;0] on e2
    Cx = Cx1+Cx2;
    Cxdata = (alpha*Cx + Vb*Ca(t)).';
    
    % Relative error
    [relerrprec,relerr] = deal(relerr,norm(Cdata-Cxdata)/norm(Cdata));
    
    % Count iteration's number
    nit = nit+1; % nit = number of iteration for a single choise of the initial parameters (a single 'attempt')
    iter = iter+1; % iter = numer of iteration for the whole algorithm (for all the 'attempt')
    
    % Store the relative error 'relerr' and the solution 'x' for each iteration 'iter'
    vect_relerr(iter) = relerr;
    cell_K{iter} = x;
    
    % If the algorithm doesn't reach the tolerance 'toll' asked:
    % - if it's the first 'attempt' than take the other choise for the
    %   initial points and restart the Gauss-Newton's algorithm.
    % - if it's the second, and last, 'attempt' than take the
    %   kinetic parameters leading to the smallest relerr found in
    %   all the iterations 'iter'.
    if  ( relerr>=1 && abs(relerr-relerrprec)<1e-3 ) || ( nit>=15 && relerr>=1 ) || ( nit>=15 && abs(relerr-relerrprec)<1e-2 ) || (nit>=nit_max)
        
        if attempt==1
            
            % Initial guess
            attempt = 2;
            nit = 0;
            [k1x,k2x,k3x,k4x,k5x,k6x,k7x] = deal(2e-1*rand(1),1e-2*rand(1),5e-2*rand(1),1*rand(1),1*rand(1),1e-2*rand(1),1e-4*rand(1)); % K6=10^2*K7
            x = [k1x;k2x;k3x;k4x;k5x;k6x;k7x];          
            
            % Solve DIRECT PROBLEM (with initial parameters)
            Ax = [[-(k3x+k4x);k4x;0],[k5x;-(k5x+k6x);k6x],[0;0;-k7x]];            
            Cx1 = concentration(k1x,Ax,Ca,[1;0;0],t0,C0,t,glnodes,glweights); %[Ca;0;0] on e1
            Cx2 = concentration(k2x,Ax,Ca,[0;1;0],t0,C0,t,glnodes,glweights); %[0;Ca;0] on e2
            Cx = Cx1+Cx2;
            Cxdata = (alpha*Cx + Vb*Ca(t)).';
            
            % Relative error
            relerr = norm(Cdata-Cxdata)/norm(Cdata);
                      
            continue
            
        else
        
            vect_relerr(vect_relerr==0) = [];
            x = cell_K{vect_relerr==min(vect_relerr)};
            k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4);
            k5x = x(5);k6x = x(6);k7x = x(7);
            
            % Solve DIRECT PROBLEM (with the parameters leading to the smallest relerr)
            Ax = [[-(k3x+k4x);k4x;0],[k5x;-(k5x+k6x);k6x],[0;0;-k7x]];          
            Cx1 = concentration(k1x,Ax,Ca,[1;0;0],t0,C0,t,glnodes,glweights); %[Ca;0;0] on e1
            Cx2 = concentration(k2x,Ax,Ca,[0;1;0],t0,C0,t,glnodes,glweights); %[0;Ca;0] on e2
            Cx = Cx1+Cx2;
            Cxdata = (alpha*Cx + Vb*Ca(t)).';
      
            relerr = min(vect_relerr);
            
            break   
            
        end
        
    end
    
%     figure(1)
%     plot(t,Cdata,'b')
%     hold on
%     plot(t,Cxdata,'c--')
% 
%     pause
      
end

