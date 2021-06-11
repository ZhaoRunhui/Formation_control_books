%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines the desired velocity
% Input variable: time t
% Output variable: desired velocity vd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vd] = Desired_2D_RT_velocity(t,q)
w = [0 0 1];
vt = [1;cos(t);0];
n = 6;
vr = zeros(n,3);
q = q';
qr = q - kron(ones(6,1),q(n,:));

for i = 1:n
    vr(i,:) = cross(w,qr(i,:));
end
% vr
vd = kron(ones(n,1),vt)+reshape(vr',[],1); 
% vd = [vd_x ; vd_y ; vd_z]