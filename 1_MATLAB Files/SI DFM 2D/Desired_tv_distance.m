%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines the desired, time-varying distances and their
% derivative
% Input variables: time t, number of agents n, and Adjacency matrix Adj
% Output variables: desired distance vector d(t) and derivative d_dot(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,d_dot] = Desired_tv_distance(t,n,Adj)

% Calculation of d(t)
A = 0.5;                    % Amplitude of AC component of desired distance
w = 0.4;                    % Frequency of AC component of desired distance 
lb = 0.5;                   % lower bound of circumradius of pentegon
AC = A*sin(w*t)+A+lb;       % AC component of desired distances
x1 = 0; % + AC;             % x coordinate of vertex 1*
y1 = AC;                    % y coordinate of vertex 1*
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
% x coordinate of framework F*(t)
x_coor = [x1; -s1*AC; -s2*AC; s2*AC; s1*AC];
% y coordinate of framework F*(t) 
y_coor = [y1; c1*AC; -c2*AC; -c2*AC; c1*AC];
q_star2 = [x_coor'; y_coor'];        % 2xn vector

d = zeros(n,n);                      % initialize the desired distance
d_dot = zeros(2*n-3,1);

for i = 1:n
    for j = 1:n
        d(i,j) = sqrt((x_coor(i)-x_coor(j))^2 + (y_coor(i)-y_coor(j))^2);
    end
end

% Calculation of d_dot(t)
ACdot = A*w*cos(w*t); 
x1dot = 0; 
% xdot coordinate of framework F*(t)
x_coor_dot = [x1dot; -s1*ACdot; -s2*ACdot; s2*ACdot; s1*ACdot];
% ydot coordinate of framework F*(t)
y_coor_dot = [ACdot; c1*ACdot; -c2*ACdot; -c2*ACdot; c1*ACdot];     
q_dot_star2 = [x_coor_dot'; y_coor_dot'];  % 2xn vector

for i = 1:n-1
    for j = i+1:n
        if Adj(i,j) == 1
            d_dot(i,j) = (q_star2(:,i)-q_star2(:,j))'...
                          *(q_dot_star2(:,i)-q_dot_star2(:,j))/d(i,j);
        end
    end
end