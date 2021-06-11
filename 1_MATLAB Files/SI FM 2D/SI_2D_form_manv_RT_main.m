%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agent formation problem
% Model Type: Single integrator 2D
% Control Objective: Formation maneuvering with rotation and tTranslation
% Number of Agents: 5 + 1 assist point 
% Subfunction: SI_2D_form_manv_RT_func and Desired_2D_RT_velocity
% Desired Formation: Pentagon with leader at geometric center
% Initial Conditions: Small pertubation around the desired formation
% You may use "SI_2D_form_manv_RT_fig_plot.m" to plot figures
% You may use "SI_2D_form_manv_RT_generate_animation.m" to generate gif 
% animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 6;                              % Number of agents
kv = 1;                           	% Control gain kv
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0 1 0 0 0 1;
       1 0 1 0 0 1;
       0 1 0 1 0 1;
       0 0 1 0 1 1;
       0 0 0 1 0 1;
       1 1 1 1 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting desired formation (A pentagon for example)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(n,n);                     % initialize the desired distance
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
x_coor = [0; -s1; -s2; s2; s1; 0];  % x coordinate of framework F*
y_coor = [1; c1; -c2; -c2; c1; 0];  % y coordinate of framework F* 

z_coor = [0; 0; 0; 0; 0; 0];        % z coordinate of framework F*
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
tfinal = 10;                        % Simulation ending time assume always
                                    % starts at 0
h = 1e-2;                           % ODE step
% Encapusulate the paremeters into a structure para
para.n = n;                         
para.kv = kv;
para.Adj = Adj;
para.d = d;
%% Randomly setting the initial conditions
% Randomly generate the initial position by a small pertubation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% q_0 = [x_coor y_coor z_coor]'+(lb+(ub-lb)*rand(3,n));

% An example initial conditions: q_0 is initial positions and q_0(:,i) is
% the coordinates for the ith agent;
q_0 = [0.3147   -0.5377   -0.8093    1.0527    1.4082   -0.3581;
       1.4058    0.4414   -0.7621   -1.1514    0.2944   -0.0782;
            0         0         0         0         0         0];

%% ODE
q_0_vec = reshape(q_0,1,[]);      	% reshape q_0 into a vector
time_span = 0:h:tfinal;             % simulation time span
[t,q] = ode45(@SI_2D_form_manv_RT_func, time_span, q_0_vec,[],para);  % ODE
% Extract information from q.
% xx(:,i) is the x coordinate for the ith agent from time 0 to tfinal;
% yy(:,i) is the y coordinate for the ith agent from time 0 to tfinal;
% zz(:,i) is the z coordinate for the ith agent from time 0 to tfinal;
xx = q(:,3*(0:n-1)+1);
yy = q(:,3*(0:n-1)+2);
zz = q(:,3*(0:n-1)+3);
%% Retrieve the control input 
u = zeros(3*n,length(t));  
% Control Input
for ii = 1:length(t)                % loop for time from 0 to tfinal
    e = zeros(n,n);                 % initialize distance error
    z = zeros(2*n-3,1);             % initialize z       
    R = zeros(2*n-3,2*n);           % initialize Rigidity Matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reconstruct R, obtain e and z 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                ord = ord+1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vd = Desired_2D_RT_velocity(t(ii),...
         reshape(q(ii,:),3,[]));  % Set desired velocity by function                                                              
    ua = -kv*R'*z;                % acquisition control input for
                                  % single integrator model
    u(:,ii) = ua+vd;              % Control input for formation manuevering
end

%% Calculate distance errors (e12,e13...)
for i = 1:n
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
            '=sqrt((xx(:,i)-xx(:,j)).^2+(yy(:,i)-yy(:,j)).^2+(zz(:,i)-zz(:,j)).^2)-d(i,j);'])
    end
end
save('SI_2D_form_manv_RT_results.mat') % save variables for plot/animation
display('Simulation complete!')