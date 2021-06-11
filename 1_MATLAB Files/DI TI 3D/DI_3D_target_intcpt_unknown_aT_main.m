%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agent Target Tracking Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Type:
%   Double integrator 3D Model
% Control Objective: 
%   Target Interception with unknown target acceleration
% Velocity Type:
%   3D Transitional velocity
% Number of Agents: 
%   8 followers + 1 leader
% Subfunction: 
%   DI_3D_target_intcpt_unknown_aT_func.m
%   Target_3D_velocity.m
% Desired Formation: 
%   Cube with the edge length being 2
% Initial Conditions: 
%   Small pertubation around the desired formation and random velocity
% You may use these files to plot figures:
%   DI_3D_target_intcpt_unknown_aT_fig_plot.m
%   DI_3D_target_intcpt_unknown_aT_results.mat
% Generate gif animation:
%   DI_3D_target_intcpt_unknown_aT_generate_animation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 9;                              % Number of agents
% Control gains kv, ka, kT and ks
kv = 0.7; 
ka = 0.2;
kT = 0.5;
ks = 3;
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0  0  0  1  1  0  0  0  1;                    
       0  0  1  0  0  1  0  0  1;
       0  1  0  1  0  0  1  0  1;
       1  0  1  0  1  0  0  1  1;
       1  0  0  1  0  1  0  1  1;
       0  1  0  0  1  0  1  1  1;
       0  0  1  0  0  1  0  1  1;
       0  0  0  1  1  1  1  0  1;
       1  1  1  1  1  1  1  1  0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting desired formation (cube with edge length = 2 for example)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(n,n);                           % initialize the desired distance
x_coor = [-1; 1; 1;-1;-1; 1; 1;-1; 0];    % x coordinate of framework F*(t)
y_coor = [-1;-1; 1; 1;-1;-1; 1; 1; 0];    % y coordinate of framework F*(t) 
z_coor = [-1;-1;-1;-1; 1; 1; 1; 1; 0];    % z coordinate of framework F*(t)
for ii = 1:n
    for jj = 1:n
        d(ii,jj) = sqrt((x_coor(ii)-x_coor(jj))^2+...
                        (y_coor(ii)-y_coor(jj))^2+...
                        (z_coor(ii)-z_coor(jj))^2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = 0.3;                           % Upper bound for random ini. condition
lb = -0.3;                          % Lower bound for random ini. condition
tfinal = 30;                        % Simulation ending time assume always
                                    % starts at 0
h = 1e-2;                           % ODE step
% Encapusulate the paremeters into a structure para
para.n = n;                         
para.kv = kv;
para.ka = ka;
para.kT = kT;
para.ks = ks;
para.Adj = Adj;
para.d = d;
%% Randomly setting the initial conditions
% Randomly generate the initial position by a small pertubation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% q_0 = [x_coor y_coor z_coor]'+(lb+(ub-lb)*rand(3,n));
% v_0 = zeros(3,n)+(lb+(ub-lb)*rand(3,n))*2;    % Initial velocity
% qv_0 = [q_0 v_0];

% An example initial conditions: q_0 is initial positions and q_0(i,:) is
% the coordinates for the ith agent; v_0 is inital velocities and v_0(i,:)
% is velocity for the ith agent; 
%%%%% initial 3D positions for n agents 3*n %%%%%
q_0 = [-1.0647    1.1236    0.7277   -0.8831   -1.2793    1.1593    0.9939   -0.8744    0.1078;
       -0.9067   -1.2809    0.7583    0.8903   -1.0368   -0.8229    0.9674    1.1528    0.0931;
       -1.1973   -1.1338   -0.8059   -0.7299    0.9289    0.8121    1.0878    0.8656   -0.2024];
%%%%% initial 3D velocities for n agents 3*n %%%%%
v_0 = [-0.4572   -0.1915    0.3015    0.2389    0.0567   -0.2910    0.3771   -0.1800    0.1393
       -0.0020    0.1023   -0.2939    0.4691   -0.4337    0.4089   -0.3078   -0.3641   -0.0321
        0.5517   -0.3314    0.0071    0.5511   -0.4208   -0.2949    0.5151   -0.2987   -0.1780];
qt_0 = [0;2;0];       % initial condition for target
qvt_0 = [q_0 v_0 qt_0];

%% ODE
qvt_0_vec = reshape(qvt_0,1,[]);      % reshape qvt_0 into a vector
time_span = 0:h:tfinal;               % simulation time span
[t,qvt] = ode45(@DI_3D_target_intcpt_unknown_aT_func,...
               time_span, qvt_0_vec,[],para);% ODE
% Extract information from qvt.
% xx(:,i) is the x coordinate for the ith agent from time 0 to tfinal;
% yy(:,i) is the y coordinate for the ith agent from time 0 to tfinal;
% zz(:,i) is the z coordinate for the ith agent from time 0 to tfinal;
% vx(:,i) is the x velocity for the ith agent from time 0 to tfinal;
% vy(:,i) is the y velocity for the ith agent from time 0 to tfinal;
% vz(:,i) is the z velocity for the ith agent from time 0 to tfinal;
% xxt is the x coordinate for the target from time 0 to tfinal;
% yyt is the y coordinate for the target from time 0 to tfinal;
% zzt is the z coordinate for the target from time 0 to tfinal;
xx = qvt(:,3*(0:n-1)+1);
yy = qvt(:,3*(0:n-1)+2);
zz = qvt(:,3*(0:n-1)+3);
vx = qvt(:,3*(0:n-1)+3*n+1);
vy = qvt(:,3*(0:n-1)+3*n+2);
vz = qvt(:,3*(0:n-1)+3*n+3);
xxt = qvt(:,6*n+1);
yyt = qvt(:,6*n+2);
zzt = qvt(:,6*n+3);

%% Retrieve the control input 
vf = zeros(3*n,length(t));          % fictitious velocity
u = zeros(3*n,length(t));           % Control Input
s = zeros(3*n,length(t));
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
            e(i,j)=sqrt((xx(ii,i)-xx(ii,j))^2+...
                   (yy(ii,i)-yy(ii,j))^2+(zz(ii,i)-zz(ii,j))^2)-d(i,j);
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eT = [xxt(ii);yyt(ii);zzt(ii)]-[xx(ii,n);yy(ii,n);zz(ii,n)];
    eTdot = Target_3D_velocity(t(ii))-[vx(ii,n);vy(ii,n);vz(ii,n)];
    v = reshape([vx(ii,:);vy(ii,:);vz(ii,:)],[],1);  % obtain v
    ua = -kv*R'*z;                      % acquisition control input for
                                        % single integrator model
    uadot = -kv*R_dot'*z-kv*R'*2*R*v;   % time derivative of ua
    vf(:,ii) = ua+kron(ones(n,1),Target_3D_velocity(t(ii))+kT*eT);   
    s(:,ii) = v-vf(:,ii);                         
    hd = uadot-ks*sign(s(:,ii))+kron(ones(n,1),kT*eTdot);
    u(:,ii) = -ka*s(:,ii)+hd-R'*z;   % u(:,ii) is the control input for all 
                                     % the agents at time ii
end

%% Calculate distance errors (e12,e13...)
for i = 1:n-1
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
            '=sqrt((xx(:,i)-xx(:,j)).^2+(yy(:,i)-yy(:,j)).^2+(zz(:,i)-zz(:,j)).^2)-d(i,j);'])
    end
end
% save variables for plot and animation
save('DI_3D_target_intcpt_unknown_aT_results.mat') 
display('Simulation complete!')