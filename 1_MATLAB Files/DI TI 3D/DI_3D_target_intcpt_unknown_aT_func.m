%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for target interception of double integrator model
% input variables: t as time sequence; qvt_vec as vector of q, v and qt;
%                  para is the structured parameters passing from the main
% Output variable: dqvt is the time derivative of q, v and qt
% A subfunction "Target_3D_velocity" is used to set the target's velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dqvt ] = DI_3D_target_intcpt_unknown_aT_func( t,qvt_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
ka = para.ka;
ks = para.ks;
kT = para.kT;
Adj = para.Adj;
d = para.d;
%
% Obtain q and v from vector qv_vec such that
% q(:,i) indicates the coordinate of the ith agent
% v(:,i) indicates the x and y direcation velocity of the ith aget
qvt = reshape(qvt_vec,3,[]);
q = qvt(:,1:n);
v = qvt(:,n+1:2*n);
qt = qvt(:,2*n+1);
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
vt = Target_3D_velocity(t); % set target's velocity
eT = qt-q(:,n);             % distance error between target and leader
eTdot = vt-v(:,n);          % time derivative of eT
v = qvt_vec(3*n+1:6*n);     % Change v from  2Xn matrix to a column vector
ua = -kv*R'*z;           % Acquisition control for single integrator model
uadot = -kv*R_dot'*z-kv*R'*2*R*v;             % time derivative for ua
vf = ua+kron(ones(n,1),vt+kT*eT);             % fictitious velocity
s = v-vf;
hd = uadot-ks*sign(s)+kron(ones(n,1),kT*eTdot);
u = -ka*s+hd-R'*z;           % control input
dqvt = [v;u;vt];