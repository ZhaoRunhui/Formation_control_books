%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for fomration acquisition of vehicle dynamic model with
% known dynamics.
% input variables: t as time sequence; q_vec as a vector of pc, eta and
%                  phihat
%                  para is the structured parameters passing from the main
% Output variable: dq is the time derivative of pc, eta and phihat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dq ] = VD_form_acq_adaptive_func( t,q_vec,para )
% Obtain parameters from structure para
n = para.n;
kv = para.kv;
ka = para.ka;
Gamma = para.Gamma;
Adj = para.Adj;
d = para.d;
L = para.L;
Mbar = para.Mbar;
Dbar = para.Dbar;
%
% Structure for q_stack
% [  xc1     xc2     ...     xcn
%    yc1     yc2     ...     ycn
%  theta1  theta2    ...   thetan
%  vcnorm1 vcnorm2   ...   vcnormn
%  omega1  omega2    ...   omegan]
q_stack = reshape(q_vec(1:5*n),5,[]);                      
pc = q_stack(1:3,:);
eta = q_stack(4:5,:);
qc = q_stack(1:2,:);
theta = (q_stack(3,:))';
vcnorm = (q_stack(4,:))';
dtheta = (q_stack(5,:))';
q = qc+kron(L',ones(2,1)).*[cos(theta) sin(theta)]';
phihat = q_vec(5*n+1:11*n);
% Find vh through theta, vcnorm and omega
dq = zeros(2,n);
M = zeros(2*n,2*n);
C = zeros(2*n,2*n);
D = zeros(2*n,2*n);
for ii = 1:n
    dq(:,ii) = [cos(theta(ii)) -L(ii)*sin(theta(ii));
                sin(theta(ii)) L(ii)*cos(theta(ii))]...
                *[vcnorm(ii);dtheta(ii)];
    J = [cos(theta(ii)) sin(theta(ii));
        -sin(theta(ii))/L(ii) cos(theta(ii))/L(ii)];
    J_dot = [-sin(theta(ii)) cos(theta(ii));
        -cos(theta(ii))/L(ii) -sin(theta(ii))/L(ii)]*dtheta(ii);
    M(2*ii-1:2*ii,2*ii-1:2*ii) = J'*Mbar(:,2*ii-1:2*ii)*J;
    C(2*ii-1:2*ii,2*ii-1:2*ii) = J'*Mbar(:,2*ii-1:2*ii)*J_dot;
    D(2*ii-1:2*ii,2*ii-1:2*ii) = J'*Dbar(:,2*ii-1:2*ii)*J;
end
%
z = zeros(2*n-3,1);      % initialize z
R = zeros(2*n-3,2*n);    % initialize Rigidity Matrix
R_dot = zeros(2*n-3,2*n);% initialize Rigidity Matrix with v_ij as its term
e = zeros(n,n);          % Distance error matrix
ord = 1;
for i = 1:n-1
    for j = i+1:n
        e(i,j) = sqrt((q(:,i)-q(:,j))'*(q(:,i)-q(:,j)))-d(i,j);
        if Adj(i,j) == 1
            z(ord) = e(i,j)*(e(i,j)+2*d(i,j));
            R(ord,2*i-1:2*i) = (q(:,i)-q(:,j))';
            R(ord,2*j-1:2*j) = (q(:,j)-q(:,i))';
            R_dot(ord,2*i-1:2*i) = (dq(:,i)-dq(:,j))';
            R_dot(ord,2*j-1:2*j) = (dq(:,j)-dq(:,i))';
            ord = ord+1;
        end
    end
end

dq = reshape(dq,[],1);     % Head velocity
vf = -kv*R'*z;             % acquisition control of single integrator model
vfdot = -kv*R_dot'*z-kv*R'*2*R*dq; % derivative of vf
s = dq-vf;
Y = regMatrix(n,theta,dtheta,dq,vf,vfdot);  % Get regression matrix
u = -ka*s+Y*phihat-R'*z;                    % Cotrol input

dpc = zeros(3,n);                           % derivative of pc
deta = zeros(2,n);                          % derivative of eta
for ii = 1:n
    S = [cos(theta(ii)) 0;sin(theta(ii)) 0;0 1];
    J = [cos(theta(ii)) sin(theta(ii));
        -sin(theta(ii))/L(ii) cos(theta(ii))/L(ii)];
    dpc(:,ii) = S*eta(:,ii);
    deta(:,ii) = inv(Mbar(:,2*ii-1:2*ii))*(inv(J')*u(2*ii-1:2*ii)-...
                 Dbar(:,2*ii-1:2*ii)*eta(:,ii));
end
dphihat = -Gamma*Y'*s;                      % estimated dynamics
dq = [reshape([dpc;deta],[],1);dphihat];