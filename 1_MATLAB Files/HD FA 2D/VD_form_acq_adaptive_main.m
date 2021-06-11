%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agent formation problem
% Model Type: Vehicle dynamic model 
% Control Objective: Formation acquisition using adaptive control
% Number of Agents: 5
% Potential Function: Quadratic 
% Subfunction: VD_form_acq_adaptive_func.m, regMatrix.m
% Desired Formation: Pentagon with the edge length being 1
% Initial Conditions: Small pertubation around the desired formation
%                     and zero velocity
% You may use "VD_form_acq_adaptive_fig_plot.m" to plot figures
% You may use "VD_form_acq_adaptive_generate_animation.m" to generate gif 
% animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 5;                            % Number of agents
kv = 1;                           % Control gain kv 
ka = 2;                           % Control gain ka
Gamma = 1*eye(6*n,6*n);           % Adaptive control gain capital gamma
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0 1 1 1 1;
       1 0 1 0 0;
       1 1 0 1 0;
       1 0 1 0 1;
       1 0 0 1 0];                
L = 0.15*ones(n,1);               % length from the mass center to hand
m = 3.6*ones(n,1);                % mass of a robot (kg)
I = 0.0405*ones(n,1);             % moment of inertia (kg-m^2)
% The equation of motion is: Mbar*eta_dot+Dbar*eta = ubar
Mbar = zeros(2,2*n);
Dbar = zeros(2,2*n);
for i = 1:n
    Mbar(:,2*i-1:2*i) = diag([m(i),I(i)]);        % mass matrix
    Dbar(:,2*i-1:2*i) = diag([0.3,0.004]);        % constant damping matrix
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting desired formation (A pentagon for example)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(n,n);                    % initialize the desired distance
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
x_coor = [0; -s1; -s2; s2; s1];    % x coordinate of framework F*
y_coor = [1; c1; -c2; -c2; c1];    % y coordinate of framework F* 
for ii = 1:n
    for jj = 1:n
        d(ii,jj) = sqrt((x_coor(ii)-x_coor(jj))^2+...
                        (y_coor(ii)-y_coor(jj))^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = 0.5;                           % Upper bound for random ini. condition
lb = -0.5;                          % Lower bound for random ini. condition
tfinal =10;                         % Simulation ending time assume always
                                    % starts at 0
h = 1e-2;                           % ODE step
% Encapusulate the paremeters into a structure para
para.n = n;                         
para.kv = kv;
para.ka = ka;
para.Gamma = Gamma;
para.Adj = Adj;
para.d = d;
para.L = L;
para.Mbar = Mbar;
para.Dbar = Dbar;
%% Randomly setting the initial conditions
% Randomly generate the initial hand position by a small pertubation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% qh_0 is initial positions and qh_0(i,:) is the coordinates for 
% "hand position" of the ith vehicle;
% qh_0 = [x_coor y_coor]'+(lb+(ub-lb)*rand(2,n)); 
% Generate random theta
% theta_0 = zeros(1,n)+(0+(2*pi-0)*rand(1,n));             % angle theta
% Obtain the coordinates for the centers of mass based on qh_0 and theta
% qc_0 = qh_0-kron(L',ones(2,1)).*[cos(theta_0);sin(theta_0)];
% Generate random center speed
% vcnorm_0 = zeros(1,n)+(lb+(ub-lb)*rand(1,n));
% Generate random angular velocity of the vehicle
% dtheta_0 = zeros(1,n)+(lb+(ub-lb)*rand(1,n));
% phihat_0 is the initial condition for the estimated dynamics 

% An example initial conditions: theta_0(i) is the initial angle of the ith
% vehicle; qc_0(i,:) is initial position of the "center of mass" of the ith
% vehicle; vcnorm_0(i) is the initial speed of the ith vehical; dtheta_0(i)
% is the initial angular velocity of the ith vehicle; phihat_0 is set to
% zero as the initial conditions of the estimated dynamics.
theta_0 = [1.1530 2.3153 3.9309 4.9023 0.5097];
qc_0 = [-0.3919   -0.6177   -0.5312    0.3558    0.5091;
        1.0120    0.3464   -0.6555   -0.4170    0.4226];
vcnorm_0 = [0.4294 0.2757 -0.0132 -0.0641 -0.0532];
dtheta_0 = [-0.1937 0.0085 0.0108 0.3176 0.2948];
phihat_0 = zeros(6*n,1);    
%% ODE
pc_0 = [qc_0;theta_0];           % d(pc_i)/dt = S(theta(i))*eta(i)
eta_0 = [vcnorm_0;dtheta_0];     
q_0_vec = [reshape([pc_0;eta_0],[],1);phihat_0]; % rearrange initial 
                                              % conditions to pass into ODE

time_span = 0:h:tfinal;
[t,q] = ode45(@VD_form_acq_adaptive_func, time_span, q_0_vec,[],para);

% Extract information from q.
% qcx(:,i) is the x coordinate for the ith vehicle center from time 0 to tfinal;
% qcy(:,i) is the y coordinate for the ith vehicle center from time 0 to tfinal;
% theta(:,i) is the hand angle of the ith vehicle from time 0 to tfinal;
% vcnorm(:,i) is the center speed of ith vehicle from time 0 to tfinal;
% dtheta(:,i) is the angular velocity of ith vehicle from time 0 to tfinal;
qcx = q(:,5*(0:n-1)+1);
qcy = q(:,5*(0:n-1)+2);
theta = q(:,5*(0:n-1)+3);
vcnorm = q(:,5*(0:n-1)+4);
dtheta = q(:,5*(0:n-1)+5);
phihat = q(:,5*n+1:11*n);
% Obtain hand coordinate and velocity according to
% qcx,qcy,theta,dtheta and vcnorm
qhx = qcx+kron(L',ones(length(t),1)).*cos(theta);
qhy = qcy+kron(L',ones(length(t),1)).*sin(theta);
vhx = vcnorm.*cos(theta)-kron(L',ones(length(t),1)).*dtheta.*sin(theta);
vhy = vcnorm.*sin(theta)+kron(L',ones(length(t),1)).*dtheta.*cos(theta);

%% Retrieve the control input
u = zeros(2*n,length(t));           % Control Input
ubar = zeros(2*n,length(t));        % Actual Control input
for ii = 1:length(t)                % loop for time from 0 to tfinal
    e = zeros(n,n);                 % initialize distance error
    z = zeros(2*n-3,1);             % initialize z       
    R = zeros(2*n-3,2*n);           % initialize Rigidity Matrix
    R_dot = zeros(2*n-3,2*n);       % initialize Rigidity Matrix 
                                    % with v_ij as its term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct R and R_dot, obtain e and z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ord = 1;
    for i = 1:n-1
        for j = i+1:n
            e(i,j) = sqrt((qhx(ii,i)-qhx(ii,j))^2+...
                          (qhy(ii,i)-qhy(ii,j))^2)-d(i,j);
            if Adj(i,j) == 1
                z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
                R(ord,2*i-1:2*i) = [qhx(ii,i)-qhx(ii,j) qhy(ii,i)-qhy(ii,j)];
                R(ord,2*j-1:2*j) = [qhx(ii,j)-qhx(ii,i) qhy(ii,j)-qhy(ii,i)];
                R_dot(ord,2*i-1:2*i) = [vhx(ii,i)-vhx(ii,j) vhy(ii,i)-vhy(ii,j)];
                R_dot(ord,2*j-1:2*j) = [vhx(ii,j)-vhx(ii,i) vhy(ii,j)-vhy(ii,i)];
                ord = ord+1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = zeros(2*n,2*n); % Mass matrix of Euler-Lagrange-like dynamics model
    C = zeros(2*n,2*n); % Coriolis matrix of Euler-Lagrange-like dynamics model
    D = zeros(2*n,2*n); % Damping matrix of Euler-Lagrange-like dynamics model
    for jj = 1:n
        J = [cos(theta(ii,jj)) sin(theta(ii,jj));
            -sin(theta(ii,jj))/L(jj) cos(theta(ii,jj))/L(jj)];
        J_dot = [-sin(theta(ii,jj)) cos(theta(ii,jj));
            -cos(theta(ii,jj))/L(jj) -sin(theta(ii,jj))/L(jj)]*dtheta(ii,jj);
        M(2*jj-1:2*jj,2*jj-1:2*jj) = J'*Mbar(:,2*jj-1:2*jj)*J;
        C(2*jj-1:2*jj,2*jj-1:2*jj) = J'*Mbar(:,2*jj-1:2*jj)*J_dot;
        D(2*jj-1:2*jj,2*jj-1:2*jj) = J'*Dbar(:,2*jj-1:2*jj)*J;
    end
    vh = reshape([vhx(ii,:);vhy(ii,:)],[],1);    
    vf = -kv*R'*z;                      % acquisition control input for
                                        % single integrator model
    vfdot = -kv*R_dot'*z-kv*R'*2*R*vh;  % time derivative of vf
    s = vh-vf;
    Y = regMatrix(n,theta(ii,:),dtheta(ii,:),vh,vf,vfdot);% Regression matrix
    % u(:,ii) is the control input for all the agents at time ii
    u(:,ii) = -ka*s+Y*(phihat(ii,:))'-R'*z; 
    for jj = 1:n
        J = [cos(theta(ii,jj)) sin(theta(ii,jj));
            -sin(theta(ii,jj))/L(jj) cos(theta(ii,jj))/L(jj)];
        ubar(2*jj-1:2*jj,ii)=inv(J')*u(2*jj-1:2*jj,ii);
    end    
end

%% Calculate distance errors (e12,e13...)
for i = 1:n-1
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
           '=sqrt((qhx(:,i)-qhx(:,j)).^2+(qhy(:,i)-qhy(:,j)).^2)-d(i,j);'])
    end
end

% %% Calculate error of estimated dynamics
phi = zeros(length(t),6*n);
for i = 1:n
    phi(:,6*i-5:6*i) = kron([m(i),I(i)/L(i)^2,Dbar(1,2*i-1),...
                            Dbar(1,2*i)/L(i), Dbar(2,2*i-1)/L(i),...
                            Dbar(2,2*i)/L(i)^2],ones(length(t),1));
end
phitilde = phihat-phi;
save('VD_form_acq_adaptive_results.mat') % save variables for plot and animation
display('Simulation complete!')