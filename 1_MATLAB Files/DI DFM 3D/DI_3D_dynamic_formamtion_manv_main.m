%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agents dynamic formation maneuvering problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Type:
%   Double integrator 3D Model (quiver tri 3D)
% Control Objective: 
%   Time Variant Distance Formation Maneuvering
% Velocity Types:
%   Angular velocity and transitional velocity
% Number of Agents: 
%   8 followers + 1 leader
% Potential Function:
%   Quadratic 
% Subfunction: 
%   DI_3D_dynamic_formamtion_manv_func.m 
%   Desired_3D_tv_RT_distance.m
%   (quiver_tri_3D.m)
% Desired Formation: 
%   cube with the edge length being 2
% Initial Conditions: 
%   Initial position is a small pertubation around the desired formation; 
%   Initial velocity is random value between 2*[lb ub]
% You may use these files to plot figures:
%   DI_3D_dynamic_formamtion_manv_fig_plot.m
% Generate gif animation:
%   DI_3D_dynamic_formamtion_manv_generate_animation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 9;                        % Number of agents
kv = 1;                       % Control gain kv 
ka = 1;                       % Control gain ka
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0 0 0 1 1 0 0 0 1;
       0 0 1 0 0 1 0 0 1;
       0 1 0 1 0 0 1 0 1;
       1 0 1 0 1 0 0 1 1;
       1 0 0 1 0 1 0 1 1;
       0 1 0 0 1 0 1 1 1;
       0 0 1 0 0 1 0 1 1;
       0 0 0 1 1 1 1 1 1;
       1 1 1 1 1 1 1 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting desired formation (cube with edge length = 2 for example)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_coor = [-1;1;1;-1;-1;1;1;-1;0];       % x coordinate of framework F*(t)
y_coor = [-1;-1;1;1;-1;-1;1;1;0];       % y coordinate of framework F*(t) 
z_coor = [-1;-1;-1;-1;1;1;1;1;0];       % z coordinate of framework F*(t) 
coor_xyz = [x_coor,y_coor,z_coor];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = 0.5;                           % Upper bound for random ini. condition
lb = -0.5;                          % Lower bound for random ini. condition
tfinal =30;                         % Simulation ending time assume always
                                    % starts at 0
h = 1e-1;                           % ODE step
% Encapusulate the paremeters into a structure para
para.n = n;                         
para.kv = kv;
para.ka = ka;
para.Adj = Adj;
para.coor_xyz = coor_xyz;
%% Randomly setting the initial conditions
% Randomly generate the initial position by a small pertubation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% q_0 = [x_coor y_coor z_coor]'+(lb+(ub-lb)*rand(3,n));
% v_0 = zeros(3,n)+(lb+(ub-lb)*rand(3,n))*2;    % Initial velocity
% qv_0 = [q_0 v_0];

% An example initial conditions: q_0 is initial positions and q_0(i,:) is
% the coordinates for the ith agent; v_0 is inital velocities and v_0(i,:)
% is x and y direction velocities for the ith agent; 
%%%%% initial 3D positions for n agents 3*n %%%%%
q_0 = [-0.6853    1.4134    0.7785   -0.5351   -0.5428    0.6419    1.2922   -1.4643    0.1787;
       -0.5942   -0.8676    1.0469    0.6576   -1.0146   -1.0782    1.4595    1.3491    0.2577;
       -1.3730   -1.4025   -0.5425   -0.5294    1.3003    1.4157    1.1557    1.4340    0.2431];
%%%%% initial 3D velocities for n agents 3*n %%%%%
v_0 = [-0.2155    0.4121   -0.9077    0.3897   -0.9311    0.5310   -0.0205    0.4187    0.3594;
        0.3110   -0.9363   -0.8057   -0.3658   -0.1225    0.5904   -0.1088    0.5094    0.3102;
       -0.6576   -0.4462    0.6469    0.9004   -0.2369   -0.6263    0.2926   -0.4479   -0.6748];
qv_0 = [q_0 v_0];                    % initial condition for the ODE solver

%% ODE
qv_0_vec = reshape(qv_0,1,[]);      % reshape qv_0 into a vector
time_span = 0:h:tfinal;             % simulation time span
[t,qv] = ode45(@DI_3D_dynamic_formamtion_manv_func,...
               time_span, qv_0_vec,[],para); % ODE
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
vf = zeros(3*n,length(t));        % fictitious velocity
u = zeros(3*n,length(t));         % Control Input - acceleration level
s = zeros(3*n,length(t)); 
for ii = 1:length(t)              % loop for time from 0 to tfinal
    e = zeros(n,n);               % initialize distance error
    z = zeros(3*n-6,1);           % initialize z       
    R = zeros(3*n-6,3*n);         % initialize Rigidity Matrix
    R_dot = zeros(3*n-6,3*n);     % initialize Rigidity Matrix differential
                                                               
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% qv(ii,:) 
    [d,dv,dv_dot,vd,vd_dot] = ...
                 Desired_3D_tv_RT_distance(t(ii),n,Adj,qv(ii,:),coor_xyz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ord = 1;                        % counter
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
    % reshape v to a column vector
    v = reshape([vx(ii,:);vy(ii,:);vz(ii,:)],[],1);  
    Rp =R'/(R*R');                   % Moore-penrose pseudoinverse 
    % derivative of M-P psedoinverse matrix
    Rp_dot = R_dot'/(R*R')-R'/(R*R')*(R_dot*R'+R*R_dot')/(R*R'); 
    % dynamic formation maneuvering control for single integrator model
    vf(:,ii) = Rp*(-kv*z+dv)+vd; 
    % time derivative of vf
    vfdot = Rp_dot*(-kv*z+dv)+Rp*(dv_dot-kv*2*(R*v-dv))+vd_dot; 
    s(:,ii) = v-vf(:,ii);         
    % u(:,ii) is the acceleration level control input for all the agents 
    % at time ii
    u(:,ii) = -ka*s(:,ii)+vfdot-R'*z; 
end

%% Calculate distance errors (e12,e13...)
d = zeros(length(t),n,n);
for ii = 1:length(t)
    [d(ii,:,:),~,~,~,~] = ...
                        Desired_3D_tv_RT_distance(t(ii),n,Adj,qv,coor_xyz);
end

for i = 1:n-1
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
            '=sqrt((xx(:,i)-xx(:,j)).^2+(yy(:,i)-yy(:,j)).^2+(zz(:,i)-zz(:,j)).^2)-d(:,i,j);'])
    end
end
% save variables for plot and animation
save('DI_3D_dynamic_formamtion_manv_results.mat') 
display('Simulation complete!')