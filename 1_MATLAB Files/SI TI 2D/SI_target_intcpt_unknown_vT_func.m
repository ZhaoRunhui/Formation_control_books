%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for target interception of single integrator model with
% unknown target velocity
% input variables: t as time sequence; qqt_vec as vector of q and qt
%                  para is the structured parameters passing from the main
% Output variable: dqqt is the time derivative of q, qt and VThat
% A subfunction "Target_velocity" is used to set up the target velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dqqt ] = SI_target_intcpt_unknown_vT_func( t,qqt_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
k1 = para.k1;
k2 = para.k2;
Adj = para.Adj;
d = para.d;
%
% Obtain q from vector qqt_vec such that
% q(:,i) indicates the coordinate of the ith agent
% qt is the absolute position of the target
% VThat is the estimated velocity of the target
qqt = reshape(qqt_vec,2,[]);
q = qqt(:,1:n);
qt = qqt(:,n+1);
VThat = qqt(:,n+2);
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
            R(ord,2*i-1:2*i) = (q(:,i)-q(:,j))';
            R(ord,2*j-1:2*j) = (q(:,j)-q(:,i))';
            ord = ord+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eT = qt-q(:,n);           % distance error between target and leader
ua = -kv*R'*z;            % Acquisition control for single integrator model
hs = kron(ones(n,1),(k1+1)*eT+VThat-ua(2*n-1:2*n));
u = ua+hs;                % control input
dqt = Target_velocity(t);   
dVThat = k1*eT+k2*sign(eT);
dqqt = [u;dqt;dVThat];