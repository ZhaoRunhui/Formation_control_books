%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agent formation problem
% Model Type: Single integrator 2D
% Control Objective: Dynamic formation maneuvering
% Number of Agents: 5 
% Subfunction: SI_dynamic_fomation_manv_func
% Desired Formation: Pentagon with the edge length being 1
% Initial Conditions: Small pertubation around the desired formation
% You may use "SI_dynamic_fomation_manv_fig_plot.m" to plot figures
% You may use "SI_dynamic_fomation_manv_generate_animation.m" to generate 
% gif animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 5;                              % Number of agents
kv = 10;                             % Control gain kv 
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0 1 1 1 1;
       1 0 1 0 0;
       1 1 0 1 0;
       1 0 1 0 1;
       1 0 0 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = 0.5;                        % Upper bound for random initial condition
lb = -0.5;                       % Lower bound for random initial condition

tfinal = 15;                     % Simulation ending time
h = 1e-1; % 1e-2                 % ODE step
% Encapusulate the paremeters into a structure 'para'
para.n = n;                         
para.kv = kv;
para.Adj = Adj;
%% Randomly setting the initial conditions
% Randomly generate the initial position by a small perturbation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% q_0 = [x_coor y_coor]'+(lb+(ub-lb)*rand(2,n));

% An example initial conditions: q_0 is initial positions and q_0(:,i) is
% the coordinates for the ith agent;
q_0 = [0.3010   -0.5222   -0.8213    0.7251    1.4141;
       0.5292    0.5393   -0.9305   -0.8502    0.3558];

%% ODE
q_0_vec = reshape(q_0,1,[]);      % reshape q_0 into a 2nx1 vector
time_span = 0:h:tfinal;             % simulation time span
[t,q] = ode45(@SI_dynamic_fomation_manv_func, time_span, q_0_vec,[],para);
% Extract information from q.
% xx(:,i) is the x coordinate for the ith agent from time 0 to tfinal;
% yy(:,i) is the y coordinate for the ith agent from time 0 to tfinal;
xx = q(:,2*(0:n-1)+1);
yy = q(:,2*(0:n-1)+2);

%% Retrieve the control input
u = zeros(2*n,length(t));           % Control Input
for ii = 1:length(t)                % loop for time from 0 to tfinal
    e = zeros(n,n);                 % initialize distance error
    dv = zeros(2*n-3,1);            % initialize dv
    z = zeros(2*n-3,1);             % initialize z
    R = zeros(2*n-3,2*n);           % initialize Rigidity Matrix
    [d,d_dot] = Desired_tv_distance(t(ii),n,Adj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct R, e, and z  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ord = 1;
    for i = 1:n-1
        for j = i+1:n
            e(i,j)=sqrt((xx(ii,i)-xx(ii,j))^2+(yy(ii,i)-yy(ii,j))^2)...
                   -d(i,j);
            if Adj(i,j) == 1
                dv(ord) = d(i,j)*d_dot(i,j);
                z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
                R(ord,2*i-1:2*i) = [xx(ii,i)-xx(ii,j) yy(ii,i)-yy(ii,j)];
                R(ord,2*j-1:2*j) = [xx(ii,j)-xx(ii,i) yy(ii,j)-yy(ii,i)];
                ord = ord+1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vd = Desired_velocity(t(ii));   % Set desired velocity by function
    % maneuvering control input for single integrator model
    u(:,ii) = R'*inv(R*R')*(-kv*z + dv)+kron(ones(n,1),vd);  
end

%% Calculate distance errors (e12,e13...)
d = zeros(length(t),n,n);
for ii = 1:length(t)
    [d(ii,:,:),~] = Desired_tv_distance(t(ii),n,Adj);
end
for i = 1:n-1
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
            '=sqrt((xx(:,i)-xx(:,j)).^2+(yy(:,i)-yy(:,j)).^2)-d(:,i,j);'])
    end
end
% save variables for plot and animation
save('SI_dynamic_fomation_manv_results.mat') 
display('Simulation complete!')