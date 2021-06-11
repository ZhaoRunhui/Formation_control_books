% Obtain regression matrix Y(q,dq,miu,dmiu)
function Y = regMatrix(n,theta,dtheta,dq,miu,dmiu)
Y = zeros(2*n,6*n);
for ii = 1:n
    sc = sin(theta(ii))*cos(theta(ii));
    s2 = sin(theta(ii))^2;
    c2 = cos(theta(ii))^2;    
    Y(2*ii-1,6*ii-5) = [c2 sc]*dmiu(2*ii-1:2*ii)...
                       +dtheta(ii)*[-sc c2]*miu(2*ii-1:2*ii);  
    Y(2*ii-1,6*ii-4) = [s2 -sc]*dmiu(2*ii-1:2*ii)...
                       +dtheta(ii)*[sc s2]*miu(2*ii-1:2*ii);
    Y(2*ii-1,6*ii-3) = [c2 sc]*dq(2*ii-1:2*ii);
    Y(2*ii-1,6*ii-2) = [-sc c2]*dq(2*ii-1:2*ii);
    Y(2*ii-1,6*ii-1) = -[sc s2]*dq(2*ii-1:2*ii);
    Y(2*ii-1,6*ii) = [s2 -sc]*dq(2*ii-1:2*ii);
    Y(2*ii,6*ii-5) = [sc s2]*dmiu(2*ii-1:2*ii)...
                     +dtheta(ii)*[-s2 sc]*miu(2*ii-1:2*ii);
    Y(2*ii,6*ii-4) = [-sc c2]*dmiu(2*ii-1:2*ii)...
                     -dtheta(ii)*[c2 sc]*miu(2*ii-1:2*ii);
    Y(2*ii,6*ii-3) = [sc s2]*dq(2*ii-1:2*ii);
    Y(2*ii,6*ii-2) = [-s2 sc]*dq(2*ii-1:2*ii);
    Y(2*ii,6*ii-1) = [c2 sc]*dq(2*ii-1:2*ii);
    Y(2*ii,6*ii) = [-sc c2]*dq(2*ii-1:2*ii);
end