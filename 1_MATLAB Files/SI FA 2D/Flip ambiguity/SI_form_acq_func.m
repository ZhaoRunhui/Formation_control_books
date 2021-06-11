%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for formation acquisition of single integrator model.
% Input variables: t as time sequence; q_vec as vector of q,
%                  para is the structured parameters passing from the main
%                  function
% Output variable: dq is the time derivative of q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dq ] = SI_form_acq_func( t,q_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
Adj = para.Adj;
d = para.d;
%
% Obtain q from vector q_vec such that
% q(:,i) indicates the coordinate of the ith agent
q = reshape(q_vec,2,[]);
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
u = -kv*R'*z;        % Acquisition control for single integrator model
dq = u;