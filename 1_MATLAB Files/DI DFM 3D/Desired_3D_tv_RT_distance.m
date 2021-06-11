%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines:
% the desired time-varying distances---------dv = dij * dij_dot
% And its derivative-----------------------------dv_dot
% actual distance between agents-------------d
% the desired swarm velocity ----------------vd
% vd's derivative ---------------------------vd_dot
% the coordinate of n vertices --------------q R(3,n)
% Input variables
%   time t, number of agents n, and Adjacency matrix Adj
% Output variables: 
%   desired distance d_bars, actual distance d and desired swarm velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,dv,dv_dot,vd,vd_dot] = ...
                            Desired_3D_tv_RT_distance(t,n,Adj,q_v,coor_xyz)

% Calculation of d(t)
A = 0.5;                    % Amplitude of AC component of desired distance
w = 0.4;                    % Frequency of AC component of desired distance 
lb = 0.5;                   % lower bound of circumradius of pentegon
AC = A*sin(w*t)+A+lb;       % AC component of desired distances
% x1 = 0; % + AC;           % x coordinate of vertex 1*
% y1 = AC;                        % y coordinate of vertex 1*
% z1 = 0;
% c1 = cos(2*pi/5);
% c2 = cos(pi/5);
% s1 = sin(2*pi/5);
% s2 = sin(4*pi/5);
x_coor = AC .* coor_xyz(:,1);     % x coordinate of framework F*(t)
y_coor = AC .* coor_xyz(:,2);     % y coordinate of framework F*(t) 
z_coor = AC .* coor_xyz(:,3);     % z coordinate of framework F*(t)
q_star3 = [x_coor'; y_coor';z_coor'];        % 3xn vector\

%

d = zeros(n,n);                      % initialize the actual distance
dv = zeros(3*n-6,1);
dv_dot = zeros(3*n-6,1);

for i = 1:n
    for j = 1:n
        d(i,j) = sqrt((x_coor(i)-x_coor(j))^2+(y_coor(i)-y_coor(j))^2+...
                 (z_coor(i)-z_coor(j))^2);
    end
end

% Calculation of d_dot(t)
ACdot = A*w*cos(w*t); 
              
x_coor_dot = ACdot .* coor_xyz(:,1);   % xdot coordinate of framework F*(t)
y_coor_dot = ACdot .* coor_xyz(:,2);   % ydot coordinate of framework F*(t)
z_coor_dot = ACdot .* coor_xyz(:,3);   % zdot coordinate of framework F*(t)
q_dot_star3 = [x_coor_dot'; y_coor_dot';z_coor_dot'];          % 3xn vector

% Calculation of d_2dot(t)
AC2dot = -A*w*w*sin(w*t);
 
x_coor_2dot = AC2dot .* coor_xyz(:,1); % xdot coordinate of framework F*(t)
y_coor_2dot = AC2dot .* coor_xyz(:,2); % ydot coordinate of framework F*(t)
z_coor_2dot = AC2dot .* coor_xyz(:,3); % zdot coordinate of framework F*(t)
q_2dot_star3 = [x_coor_2dot'; y_coor_2dot';z_coor_2dot'];      % 3xn vector

index_e = 1;
for i = 1:n-1
    for j = i+1:n
        if Adj(i,j) == 1
            dv(index_e) = (q_star3(:,i)-q_star3(:,j))'*...
                          (q_dot_star3(:,i)-q_dot_star3(:,j));
            dv_dot(index_e) = (q_star3(:,i)-q_star3(:,j))'*...
                              (q_2dot_star3(:,i)-q_2dot_star3(:,j))+...
                              (q_dot_star3(:,i)-q_dot_star3(:,j))'*...
                              (q_dot_star3(:,i)-q_dot_star3(:,j));
            index_e = index_e + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qv = reshape(q_v,3,[]);
q = qv(:,1:n);
v = qv(:,n+1:2*n);
q = q';
v = v';
ww = [1 1 1]';                                   % anglular velocity 
vt  = [1,cos(t),0]';                             % translation velocity
vt_dot = [0,-sin(t),0]';                         % translation acceleration
%%% relative position and velocity with respect to leader(agent n) %%%
qr = q - kron(ones(n,1),q(n,:));
qvr = v - kron(ones(n,1),v(n,:));
%%% Angular velocity based on relative postion %%%
vr = zeros(n,3);
vr_dot = zeros(n,3);
for ii = 1:n
    vr(ii,:) = cross(ww,qr(ii,:));
    vr_dot(ii,:) = 0 + cross(ww,qvr(ii,:));       % ww is constant here
end
%%% add vt and vr up %%%
vd  = kron(ones(n,1),vt) + reshape(vr',[],1);
vd_dot = kron(ones(n,1),vt_dot) + reshape(vr_dot',[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% q_v will be updated after every loop %%%%




