%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for 3D fomration acquisition of double integrator model
% input variables: t as time sequence; qv_vec as vector of q and v combined
%                  para is the structured parameters passing from the main
% Output variable: dqv is a vector with each time to be the time derivative
%                  of q and v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dqv ] = DI_3D_form_acq_func( t,qv_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
Adj = para.Adj;
d = para.d;
ka = para.ka;
%
% Obtain q and v from vector qv_vec such that
% q(:,i) is the coordinate of the ith agent
% v(:,i) is the velocity of the ith agent
qv = reshape(qv_vec,3,[]);
q = qv(:,1:n);
v = qv(:,n+1:2*n);
%
z = zeros(3*n-6,1);             % initialize z
R = zeros(3*n-6,3*n);           % initialize Rigidity Matrix
R_dot = zeros(3*n-6,3*n);       % Rigidity Matrix with v_ij as its term
e = zeros(n,n);                 % Distance error matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct R and R_dot, obtain e and z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ord = 1;
for i = 1:n-1
    for j = i+1:n
        e(i,j) = sqrt((q(:,i)-q(:,j))'*(q(:,i)-q(:,j)))-d(i,j);
        if Adj(i,j) == 1
            z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
            R(ord,3*i-2:3*i) = (q(:,i)-q(:,j))';
            R(ord,3*j-2:3*j) = (q(:,j)-q(:,i))';
            R_dot(ord,3*i-2:3*i) = (v(:,i)-v(:,j))';
            R_dot(ord,3*j-2:3*j) = (v(:,j)-v(:,i))';
            ord = ord+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = qv_vec(3*n+1:6*n);       % Change v from  2Xn matrix to a column vector
ua = -kv*R'*z;            % Acquisition control for single integrator model
vf = ua;
vfdot = -kv*R_dot'*z-kv*R'*2*R*v; % time derivative for ua
s = v-vf;                               
u = -ka*s+vfdot-R'*z;               % control input
dqv = [v;u];