%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agents 3D formation acquisition problem
% Model Type: 3D Model & Double integrator
% Control Objective: Formation Acquisition
% Number of Agents: 8
% Subfunction: DI_3D_form_acq_func.m
% Desired Formation: Unit cube with the edge length being 2
% Initial Conditions: Initial position is a small pertubation around the
%                     desired formation; Initial velocity is random value 
%                     between 2*[lb ub]
% You may use "DI_3D_form_acq_fig_plot.m" to plot figures
% You may use "DI_3D_form_acq_generate_animation.m" to generate gif animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 8;                              % Number of agents
kv = 1;                             % Control gain kv 
ka = 1;                             % Control gain ka
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0 1 1 1 1 0 0 1;
       1 0 1 0 1 1 0 0;
       1 1 0 1 0 1 1 0;
       1 0 1 0 0 0 1 1;
       1 1 0 0 0 1 0 1;
       0 1 1 0 1 0 1 1;
       0 0 1 1 0 1 0 1;
       1 0 0 1 1 1 1 0];            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting desired formation (A unit cube for example)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(n,n);                     % initialize the desired distance
x_coor = [-1;1;1;-1;-1;1;1;-1];   % x coordinate of framework F*(t)
y_coor = [-1;-1;1;1;-1;-1;1;1];   % y coordinate of framework F*(t) 
z_coor = [-1;-1;-1;-1;1;1;1;1];   % z coordinate of framework F*(t)
for ii = 1:n
    for jj = 1:n
        d(ii,jj) = sqrt((x_coor(ii)-x_coor(jj))^2+...
                        (y_coor(ii)-y_coor(jj))^2+...
                        (z_coor(ii)-z_coor(jj))^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = 0.5;                           % Upper bound for random ini. condition
lb = -0.5;                          % Lower bound for random ini. condition
tfinal =4;                          % Simulation ending time assume always
                                    % starts at 0
h = 1e-1;                           % ODE step
% Encapusulate the paremeters into a structure para
para.n = n;                         
para.kv = kv;
para.ka = ka;
para.Adj = Adj;
para.d = d;
%% Randomly setting the initial conditions
% Randomly generate the initial position by a small pertubation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% q_0 = [x_coor y_coor z_coor]'+(lb*ones(3,n)+(ub-lb)*rand(3,n));
% v_0 = (lb*rand(3,n)+(ub-lb)*rand(3,n))*2;    % Initial velocity
% qv_0 = [q_0 v_0];

% An example initial conditions: q_0 is initial positions and q_0(i,:) is
% the coordinates for the ith agent; v_0 is inital velocities and v_0(i,:)
% is x and y direction velocities for the ith agent; qv_0 is the initial
% conditions that will pass into the ODE solver
q_0 = [-0.6853    1.4134    0.7785   -0.5351   -0.5428    0.6419    1.2922   -1.4643;
       -0.5942   -0.8676    1.0469    0.6576   -1.0146   -1.0782    1.4595    1.3491;
       -1.3730   -1.4025   -0.5425   -0.5294    1.3003    1.4157    1.1557    1.4340];
v_0 = [0.3575   -0.2155    0.4121   -0.9077    0.3897   -0.9311    0.5310   -0.0205;
       0.5155    0.3110   -0.9363   -0.8057   -0.3658   -0.1225    0.5904   -0.1088;
       0.4863   -0.6576   -0.4462    0.6469    0.9004   -0.2369   -0.6263    0.2926];
qv_0 = [q_0 v_0];

%% ODE
qv_0_vec = reshape(qv_0,1,[]);      % reshape qv_0 into a vector
time_span = 0:h:tfinal;             % simulation time span
[t,qv] = ode45(@DI_3D_form_acq_func, time_span, qv_0_vec,[],para);    % ODE
% Extract information from qv.
% xx(:,i) is the x coordinate for the ith agent from time 0 to tfinal;
% yy(:,i) is the y coordinate for the ith agent from time 0 to tfinal;
% zz(:,i) is the z coordinate for the ith agent from time 0 to tfinal;
% vx(:,i) is the x velocity for the ith agent from time 0 to tfinal;
% vy(:,i) is the y velocity for the ith agent from time 0 to tfinal;
% vz(:,i) is the z velocity for the ith agent from time 0 to tfinal;
xx = qv(:,3*(0:n-1)+1);
yy = qv(:,3*(0:n-1)+2);
zz = qv(:,3*(0:n-1)+3);
vx = qv(:,3*(0:n-1)+3*n+1);
vy = qv(:,3*(0:n-1)+3*n+2);
vz = qv(:,3*(0:n-1)+3*n+3);
%% Retrieve the control input
vf = zeros(3*n,length(t));          % fictitious velocity
u = zeros(3*n,length(t));           % Control Input
for ii = 1:length(t)                % loop for time from 0 to tfinal
    e = zeros(n,n);                 % initialize distance error
    z = zeros(3*n-6,1);             % initialize z       
    R = zeros(3*n-6,3*n);           % initialize Rigidity Matrix
    R_dot = zeros(3*n-6,3*n);       % initialize Rigidity Matrix 
                                    % with v_ij as its term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct R and R_dot, obtain e and z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ord = 1;
    for i = 1:n-1
        for j = i+1:n
            e(i,j)=sqrt((xx(ii,i)-xx(ii,j))^2+(yy(ii,i)-yy(ii,j))^2+...
                        (zz(ii,i)-zz(ii,j))^2)-d(i,j);
            if Adj(i,j) == 1
                z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
                R(ord,3*i-2:3*i) = [xx(ii,i)-xx(ii,j) yy(ii,i)-yy(ii,j)...
                                    zz(ii,i)-zz(ii,j)];
                R(ord,3*j-2:3*j) = [xx(ii,j)-xx(ii,i) yy(ii,j)-yy(ii,i)...
                                    zz(ii,j)-zz(ii,i)];
                R_dot(ord,3*i-2:3*i)=[vx(ii,i)-vx(ii,j)...
                                      vy(ii,i)-vy(ii,j)...
                                      vz(ii,i)-vz(ii,j)];
                R_dot(ord,3*j-2:3*j)=[vx(ii,j)-vx(ii,i)...
                                      vy(ii,j)-vy(ii,i)...
                                      vz(ii,j)-vz(ii,i)];
                ord = ord+1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = reshape([vx(ii,:);vy(ii,:);vz(ii,:)],[],1);  % obtain v
    ua = -kv*R'*z;                          % acquisition control input for
                                            % single integrator model
    vf(:,ii) = ua;
    vfdot = -kv*R_dot'*z-kv*R'*2*R*v;       % time derivative of ua
    s = v-vf(:,ii);                         
    u(:,ii) = -ka*s+vfdot-R'*z;     % u(:,ii) is the control input for all 
                                    % the agents at time ii
end

%% Calculate distance errors (e12,e13...)
for i = 1:n-1
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
            '=sqrt((xx(:,i)-xx(:,j)).^2+(yy(:,i)-yy(:,j)).^2+(zz(:,i)-zz(:,j)).^2)-d(i,j);'])
    end
end
save('DI_3D_form_acq_results.mat') % save variables for plot and animation
display('Simulation complete!')
h_msg = msgbox('Operation Completed', 'Success');