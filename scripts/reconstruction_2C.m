function [k1x,k2x,k3x,k4x,relerr,nit,Cx,Cxdata] = reconstruction_2C(Cdata,Ca,t,t0,C0)

% 2-COMPARTMENT C1,C2 COMPARTMENTAL SYSTEM with SINGLE INPUT FUNCTION Ca
%       K1    ______    K3    ______
%      --->  |      |  --->  |      |   
%  Ca   K2   |  C1  |   K4   |  C2  |   
%      <---  |______|  <---  |______|
%             C1=free      C2=metabolized
% (Ca=arterial IF)

% The ODE system is C'=A*C+[K1*Ca;0] with initial condition C(0)=[0;0]
% The solution is C(t)=K1*int( expm((t-u)A)*Ca(u) ,u=t0..t) [e1]

%% Initial guess
Vb = 0.2;
[k1,k2,k3,k4] = deal(5e-1,5e-1,5e-1,5e-1); 

alpha = [1-Vb,1-Vb];

% Gauss-Legendre quadrature
ngl = 4; [glnodes,glweights] = gauss_legendre(ngl);

% Initial guess
Kx = num2cell([k1,k2,k3,k4]+(rand(1,4)-4/5).*[k1,k2,k3,k4]);
[k1x,k2x,k3x,k4x] = deal(Kx{:});

x = [k1x;k2x;k3x;k4x];

% Solve DIRECT PROBLEM (with initial parameters)
Ax = [[-(k3x+k2x);k3x],[k4x;-k4x]];
Cx = concentration(k1x,Ax,Ca,[1;0],t0,C0,t,glnodes,glweights);
% Cxdata = (alpha*Cx).'; % Vb = 0
Cxdata = (alpha*Cx + Vb*Ca(t)).'; % Vb ~= 0

% Relative error
relerr = norm(Cdata-Cxdata)/norm(Cdata);

%% Regularized GAUSS - NEWTON method 

% Iteration numbers
nit = 0;
iter = 0;
cont = 0;
nit_max = 10; 

% Initializations
vect_relerr = zeros(3*nit_max,1);
cell_K = cell(1,3*nit_max);

% Lower bound for the relative error
toll = sqrt(sum(Cdata))/norm(Cdata) * 3e-1;

while relerr>toll
    
    % Count iteration's number
    nit = nit+1; % nit = number of iteration for a single choise of the initial parameters 
    iter = iter+1; % iter = numer of iteration for the whole algorithm
    
    % Derivative with respect to k1
    aux1 = dF_dk(Ax,Ca,[1;0],t0,t,glnodes,glweights,alpha).';
    % Derivative with respect to A
    aux2 = dF_dA(k1x,Ax,Ca,[1;0],t0,C0,t,glnodes,glweights,alpha).';
    D = [aux1 aux2];
    D = [D(:,1),-D(:,2),D(:,3)-D(:,2),D(:,4)-D(:,5)];
    
    % Regularization parameter --------------------------------------------
    vect_lambda = 1e4:5e2:1e6; %values for lambda to compute V_lambda
    vect_lambda = vect_lambda';    
    % Generalized Cross Validation
    [r,~] = GCV(D,Cdata-Cxdata,vect_lambda);
    %----------------------------------------------------------------------
    % '\' solve the linear system Dh=Z with Z=Cdata-Cxdata h=increment
    % ==> Tychonov regularization: (rI+D'*D)h=D'*Z
    h = (r*diag([1,1,1,1])+D.'*D)\(D.'*(Cdata-Cxdata));
    
    % Check for positive coefficients
    xph = x+h;
    if any(xph<=0)
        h(xph<=0)=0;
    end
    x = x+h;
    
    % Refresh the parameters
    k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4);
    
    % Solve DIRECT PROBLEM  --> refresh data
    Ax = [[-(k3x+k2x);k3x],[k4x;-k4x]];
    Cx = concentration(k1x,Ax,Ca,[1;0],t0,C0,t,glnodes,glweights);
    % Cxdata = (alpha*Cx).'; % Vb = 0
    Cxdata = (alpha*Cx + Vb*Ca(t)).'; % Vb ~= 0
    
    % Relative error
    % [relerrprec,relerr] = deal(relerr,norm(Cdata-Cxdata)/norm(Cdata));
    relerr = norm(Cdata-Cxdata)/norm(Cdata);
    
    % Store the relative error 'relerr' and the solution 'x' for each iteration 'iter'
    vect_relerr(iter) = relerr;
    cell_K{iter} = x;
    
    %.....................................................................%
    if  nit>=nit_max
        % ( relerr>=0.5 && abs(relerr-relerrprec)<1e-3 ) || ( nit>=10 && relerr>=0.1 ) || ( nit>=15 && abs(relerr-relerrprec)<1e-4 ) || (nit>=nit_max)
        
        nit = 0;
        cont = cont+1;
        
        if cont == 3
        
            vect_relerr(vect_relerr==0) = [];
            x = cell_K{vect_relerr==min(vect_relerr)};
            k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4);
            
            % Solve DIRECT PROBLEM (with the parameters leading to the smallest relerr)
            Ax = [[-(k3x+k2x);k3x],[k4x;-k4x]];
            Cx = concentration(k1x,Ax,Ca,[1;0],t0,C0,t,glnodes,glweights);
            % Cxdata = (alpha*Cx).'; % Vb = 0
            Cxdata = (alpha*Cx + Vb*Ca(t)).'; % Vb ~= 0
            
            relerr = min(vect_relerr);
            
            break
        else
            
            % Initial guess
            Kx = num2cell([k1,k2,k3,k4]+(rand(1,4)-4/5).*[k1,k2,k3,k4]);
            [k1x,k2x,k3x,k4x] = deal(Kx{:});
            
            x = [k1x;k2x;k3x;k4x];
            
            % Solve DIRECT PROBLEM (with initial parameters)
            Ax = [[-(k3x+k2x);k3x],[k4x;-k4x]];
            Cx = concentration(k1x,Ax,Ca,[1;0],t0,C0,t,glnodes,glweights);
            % Cxdata = (alpha*Cx).'; % Vb = 0
            Cxdata = (alpha*Cx + Vb*Ca(t)).'; % Vb ~= 0
            
            % Relative error
            relerr = norm(Cdata-Cxdata)/norm(Cdata);
            
            continue
        end
        
    end
    
%     %.....................................................................%
%     if cont==2
%         
%         cont = 0;
%         nit = 0;
%         
%         toll = toll+toll*2e-1;
%         
%         % Initial guess
%         Kx = num2cell([k1,k2,k3,k4]+(rand(1,4)-4/5).*[k1,k2,k3,k4]);
%         [k1x,k2x,k3x,k4x] = deal(Kx{:});
%         
%         x = [k1x;k2x;k3x;k4x];
%         
%         % Solve DIRECT PROBLEM (with initial parameters)
%         Ax = [[-(k3x+k2x);k3x],[k4x;-k4x]];
%         Cx = concentration(k1x,Ax,Ca,[1;0],t0,C0,t,glnodes,glweights);
%         % Cxdata = (alpha*Cx).'; % Vb = 0
%         Cxdata = (alpha*Cx + Vb*Ca(t)).'; % Vb ~= 0
%         
%         % Relative error
%         relerr = norm(Cdata-Cxdata)/norm(Cdata);
% 
%         continue
%         
%     end
    
%     figure(2)
%     plot(t,Cdata,'b')
%     hold on
%     plot(t,Cxdata,'c--')
% 
%     pause

end

end
