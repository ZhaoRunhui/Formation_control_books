%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for fomration maneuvering of single integrator model
% with rotation
% input variables: t as time sequence; q_vec as vector of q,
%                  para is the structured parameters passing from the main
% Output variable: dq is a the time derivative of q
%  "Desired_2D_RT_velocity" is used to set up the desired velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dq ] = SI_2D_form_manv_RT_func( t,q_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
Adj = para.Adj;
d = para.d;
%
% Obtain q from vector q_vec such that
% q(:,i) indicates the coordinate of the ith agent
q = reshape(q_vec,3,[]);
%
z = zeros(2*n-3,1);             % initialize z
R = zeros(2*n-3,2*n);           % initialize Rigidity Matrix
e = zeros(n,n);                 % Distance error matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct R, obtain e and z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ord = 1;
for i = 1:n-1
    for j = i+1:n
        e(i,j) = sqrt((q(:,i)-q(:,j))'*(q(:,i)-q(:,j)))-d(i,j);
        if Adj(i,j) == 1
            z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
            R(ord,3*i-2:3*i) = (q(:,i)-q(:,j))';
            R(ord,3*j-2:3*j) = (q(:,j)-q(:,i))';
            ord = ord+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vd = Desired_2D_RT_velocity(t,q);     % desired velocity
ua = -kv*R'*z;   % Acquisition control for single integrator model
u = ua+vd;   		 % Control input for formation maneuvering
dq = u;