%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-agents formation maneuvering problem ODE function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model:
%   3D,Double-integrator model; time-variant distance, Rotation &
%   translation
% input:
%   t ------------------- time sequence
%   qv_vec -------------- vector of q and v combined
%   para ---------------- the structured parameters passing from the main
% Output: 
%   dqv ----------------- a vector with each time to be the time derivative
%                         of q and v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dqv ] = DI_3D_dynamic_formamtion_manv_func( t,qv_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
Adj = para.Adj;
coor_xyz = para.coor_xyz;
ka = para.ka;
%
% Obtain q and v from vector qv_vec such that
% q(:,i) indicates the coordinate of the ith agent
% v(:,i) indicates the 3D velocity of the ith agent
qv = reshape(qv_vec,3,[]);
q = qv(:,1:n);
v = qv(:,n+1:2*n);
%
z = zeros(3*n-6,1);         % initialize z
R = zeros(3*n-6,3*n);           % initialize Rigidity Matrix (q_ij)
R_dot = zeros(3*n-6,3*n);       % Rigidity Matrix with v_ij as its term
e = zeros(n,n);                 % Distance error matrix
%
q_v = reshape(qv,1,[]);
[d,dv,dv_dot,vd,vd_dot] = Desired_3D_tv_RT_distance(t,n,Adj,q_v,coor_xyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
v = reshape(v,[],1);               % reshape v to a column vector
Rp =R'/(R*R');                     % Moore-penrose pseudoinverse 
% dynamic formation maneuvering control for single integrator model
vf =Rp*(-kv*z+dv)+kron(ones(1,1),vd);  
% derivative on Moore-penrose pseudoinverse
Rp_dot = R_dot'/(R*R')-R'/(R*R')*(R_dot*R'+R*R_dot')/(R*R'); 
% time derivative of vf
vfdot = Rp_dot*(-kv*z+dv)+Rp*(dv_dot-kv*2*(R*v-dv))+kron(ones(1,1),vd_dot); 
s = v-vf;                         
u = -ka*s+vfdot-R'*z; % control input
dqv = [v;u];