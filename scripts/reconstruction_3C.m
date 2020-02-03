function [k1x,k2x,k3x,k4x,k5x,k6x,k7x,relerr,nit,Cxdata] = reconstruction_3C(Cdata,Ca,t)
global analysis_gui

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
% The solution is C(t)=K1*int( expm((t-u)A)*Ca(u) , u=t0..t) [e1] + ...
%                        K2*int( expm((t-u)A)*Ca(u) , u=t0..t) [e2]

%% Inizial guessed parameters for Gauss-Newton's algorithm

% Fraction of blood volume
Vb = str2num(get(analysis_gui.Vb,'string'));
alpha = [1-Vb,1-Vb,1-Vb]; % --> need to sum different compartment C_.. and to obtain C_tot

% Initial values
t0 = 0; C0 = [0;0;0];

% Initial guess
[k1x,k2x,k3x,k4x,k5x,k6x] = deal(1+3*rand(1),2*rand(1),4*rand(1),1,1,0.1); 
[k1x,k2x,k3x,k4x,k5x,k6x,k7x] = deal(5e-3,5e-3,5e-2,3,3,1e-3,1e-5); % K6=10^2*K7

k7x = 10^(-2)*k6x; % physiological information
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
cont = 0;
nit_max = 50;

% Lower bound for the relative error
toll = sqrt(sum(Cdata))/norm(Cdata);

while relerr>toll || nit==0 % <-- I want the algorithm to make at least one iteration! 
    
    % Count iteration's number
    nit = nit+1; % nit = number of iteration for a single choise of the initial parameters (a single 'attempt')
    iter = iter+1; % iter = numer of iteration for the whole algorithm (for all the 'attempt')
    
    % Derivative with respect to k1
    aux1 = dF_dk(Ax,Ca,[1;0;0],t0,t,glnodes,glweights,alpha).';
    % Derivative with respect to k2
    aux2 = dF_dk(Ax,Ca,[0;1;0],t0,t,glnodes,glweights,alpha).';
    % Derivative with respect to A
    aux3 = dF_dA(k1x,Ax,Ca,[1;0;0],t0,C0,t,glnodes,glweights,alpha).';
    aux4 = dF_dA(k2x,Ax,Ca,[0;1;0],t0,C0,t,glnodes,glweights,alpha).';
    D = [aux1 aux2 aux3+aux4];
    D = [D(:,1),D(:,2),-D(:,3),D(:,4)-D(:,3),D(:,6)-D(:,7),D(:,8)-D(:,7),-D(:,11)];
    % D=[dalphaCdK1,dalphaCdK2,dalphaCdK3,dalphaCdK4,dalphaCdK5,dalphaCdK6,dalphaCdK7]
    % considering Ax=[[-(K3x+K4x);K4x;0],[K5x;-(K5x+K6x);K6x],[0;0;-K7x]];
        
    % Regularization parameter --------------------------------------------   
    vect_lambda = 1e3:5e2:1e5; %values for lambda to compute V_lambda
    vect_lambda = vect_lambda';    
    % Generalized Cross Validation
    [r,~] = GCV(D,Cdata-Cxdata,vect_lambda); 
    %----------------------------------------------------------------------
     % '\' solve the linear system Mh=D with D=Cdata-Cxdata h=increment
    % ==> Tychonov regularization: (rI+M'*M)h=M'*D
    h = (r*diag([1,1,1,1,1,1,1])+D.'*D)\(D.'*(Cdata-Cxdata));
    
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
    
    % If the algorithm doesn't reach the tolerance 'toll' asked:
    % - if it's the first 'attempt' than take the other choise for the
    %   initial points and restart the Gauss-Newton's algorithm.
    % - if it's the second, and last, 'attempt' than take the
    %   kinetic parameters leading to the smallest relerr found in
    %   all the iterations 'iter'.
    if  ( relerr>=1 && abs(relerr-relerrprec)<1e-3 ) || ( nit>=10 && relerr>=1 ) || ( nit>=25 && abs(relerr-relerrprec)<1e-3 ) || (nit>=nit_max)
                
        cont = cont+1;
        nit = 0;
        
        % Initial guess
        [k1x,k2x,k3x,k4x,k5x,k6x] = deal(1+3*rand(1),2*rand(1),4*rand(1),1,1,0.1);
        k7x = 10^(-2)*k6x; % physiological information
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
        
    end
    
    %.....................................................................%
    if cont>=5
        
        cont = 0;
        nit = 0;
        
        toll = toll+toll*2e-1;
        
        [k1x,k2x,k3x,k4x,k5x,k6x] = deal(1+3*rand(1),2*rand(1),4*rand(1),1,1,0.1);
        k7x = 10^(-2)*k6x; % physiological information
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
        
    end
      
end
