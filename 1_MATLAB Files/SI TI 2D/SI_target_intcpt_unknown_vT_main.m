%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agent formationm problem
% Model Type: Single integrator 2D
% Control Objective: Target interception with unknown target velocity
% Number of Agents: 6 (5 followers and 1 leader)
% Subfunction: SI_target_intcpt_unknown_vT_func and Target_velocity
% Desired Formation: Pentagon with the edge length being 1, the sixth agent
% in the geometric center of the pentagon
% Initial Conditions: Small pertubation around the desired formation
% You may use "SI_target_intcpt_unknown_vT_fig_plot.m" to plot figures
% You may use "SI_target_intcpt_unknown_vT_generate_animation.m" to 
% generate gif animation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Parameter Setting 
n = 6;                              % Number of agents
% Control gains kv, k1 and k2
kv = 1;                           	 
k1 = 2;
k2 = 2;
% n x n adjacency matrix for the graph. Adj(i,j) == 0, if ith agent and jth
% agent are not connected; Adj(i,j) == 1 if ith and jth agent are connected
% Change Adj if the connecting topology is changed
Adj = [0 1 0 0 1 1;
       1 0 1 0 0 1;
       0 1 0 1 0 1;
       0 0 1 0 1 1;
       1 0 0 1 0 0;
       1 1 1 1 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting desired formation (A pentagon for example)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = zeros(n,n);                     % initialize the desired distance
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
x_coor = [0; -s1; -s2; s2; s1; 0];     % x coordinate of framework F*
y_coor = [1; c1; -c2; -c2; c1; 0];     % y coordinate of framework F* 
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
h = 1e-1; %1e-2                     % ODE step
% Encapusulate the paremeters into a structure para
para.n = n;                         
para.kv = kv;
para.k1 = k1;
para.k2 = k2;
para.Adj = Adj;
para.d = d;
%% Randomly setting the initial conditions
% Randomly generate the initial position by a small pertubation
% around the desired shape with a upper bound "ub" and a lower bound "lb"
% q_0 = [x_coor y_coor]'+(lb+(ub-lb)*rand(2,n));


% An example initial conditions: q_0 is initial positions and q_0(i,:) is
% the coordinates for the ith agent, q_0(n,:) is the initial coordinate for
% the leader; qt_0 is the initial condition for the target. VThat_0 is the
% initial condition for the target's estimated velocity
q_0 = [0.3143 -0.5218 -0.8912  0.7038 0.8027 0.0853;
       0.7435  0.1590 -1.0580 -0.8357 0.6398 0.0497];
qt_0 = [1;-1]; %[2;0];       % initial condition for target
VThat_0 = [0;0];
qqt_0 = [q_0 qt_0 VThat_0];

%% ODE
qqt_0_vec = reshape(qqt_0,1,[]);      	% reshape qqt_0 into a vector
time_span = 0:h:tfinal;             	% simulation time span
[t,qqt] = ...
      ode45(@SI_target_intcpt_unknown_vT_func,time_span,qqt_0_vec,[],para);
% Extract information from qqt.
% xx(:,i) is the x coordinate for the ith agent from time 0 to tfinal;
% yy(:,i) is the y coordinate for the ith agent from time 0 to tfinal;
% xxt is the x coordinate for the target from time 0 to tfinal;
% yyt is the y coordinate for the target from time 0 to tfinal;
% VThatx is the estimated velocity of the target on the x direction
% VThaty is the estimated velocity of the target on the y direction
xx = qqt(:,2*(0:n-1)+1);
yy = qqt(:,2*(0:n-1)+2);
xxt = qqt(:,2*n+1);
yyt = qqt(:,2*n+2);
VThatx = qqt(:,2*n+3);
VThaty = qqt(:,2*n+4);
%% Retrieve the control input 
u = zeros(2*n,length(t));           % Control Input
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
            e(i,j)=sqrt((xx(ii,i)-xx(ii,j))^2+...
                       (yy(ii,i)-yy(ii,j))^2)-d(i,j);
            if Adj(i,j) == 1
                z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
                R(ord,2*i-1:2*i) = [xx(ii,i)-xx(ii,j) yy(ii,i)-yy(ii,j)];
                R(ord,2*j-1:2*j) = [xx(ii,j)-xx(ii,i) yy(ii,j)-yy(ii,i)];
                ord = ord+1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eT =[xxt(ii);yyt(ii)]-[xx(ii,n);yy(ii,n)];
    ua = -kv*R'*z;                      	% acquisition control input for
                                            % single integrator model
    % u(:,ii) is the control input for all the agents at time ii
    hs = kron(ones(n,1),(k1+1)*eT+[VThatx(ii);VThaty(ii)]-ua(2*n-1:2*n));
    u(:,ii) = ua+hs; 
end

%% Calculate distance errors (e12,e13...)
for i = 1:n-1
    for j = i+1:n
        eval(['e',int2str(i),int2str(j),...
            '=sqrt((xx(:,i)-xx(:,j)).^2+(yy(:,i)-yy(:,j)).^2)-d(i,j);'])
    end
end
save('SI_target_intcpt_unknown_vT_results.mat') % save variables for plot 
                                              	% and animation
display('Simulation complete!')
