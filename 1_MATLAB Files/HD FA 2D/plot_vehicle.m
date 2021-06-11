%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used for ploting a vehicle
% input variables: xc, yc as the coordinate of the "center of mass" of the
%                  vehicle; theta as the angle of the vehicle's hand; L is
%                  the length of the hand;
% Output variable: h is the handle object of the plotted vehicle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plot_vehicle(xc,yc,theta,L)
w = 1/5*L;
h = 1/3*L;
xh = xc+L*cos(theta);
yh = yc+L*sin(theta);
xpolybase = [0     -1/3*L  -1/3*L  1/3*L  1/3*L;
         2/3*L 3/10*L   -0.3*L  -0.3*L 3/10*L];

% Clockwise rotation matrix
Rot = [cos(pi/2-theta) sin(pi/2-theta);-sin(pi/2-theta) cos(pi/2-theta)]; 

xbaseA = [-1/3*L -1/3*L-w   -1/3*L-w  -1/3*L;
          h/2    h/2          -h/2     -h/2];
xbaseB = [1/3*L 1/3*L+w   1/3*L+w    1/3*L;
          h/2    h/2          -h/2     -h/2];
xA = zeros(2,4);
xB = zeros(2,4);
xpoly = zeros(2,5);
for ii = 1:4
    xA(:,ii) = Rot*xbaseA(:,ii)+[xc;yc];
    xB(:,ii) = Rot*xbaseB(:,ii)+[xc;yc];
end
for ii = 1:5
    xpoly(:,ii) = Rot*xpolybase(:,ii)+[xc;yc];
end
h1 = fill(xA(1,:),xA(2,:),'b','EdgeColor','k','LineWidth',2);
h2 = fill(xB(1,:),xB(2,:),'b','EdgeColor','k','LineWidth',2);
h3 = patch(xpoly(1,:),xpoly(2,:),'r');
set(h3,'FaceColor','none','LineWidth',2)
h4 = line([xc;xh],[yc;yh],'LineWidth',4);
h = [h1;h2;h3;h4];
