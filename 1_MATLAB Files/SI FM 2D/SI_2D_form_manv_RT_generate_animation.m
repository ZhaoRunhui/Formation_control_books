%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates gif animation for single integrator model
% formation maneuvering in 2D with Rotation & Translation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_2D_form_manv_RT_results.mat'
lw = 2;     						% line width
ms = 12;    						% MakerSize
fs = 18;    						% Fontsize
pausetime = 0.1;  					% Pause time
gifName = 'SI_2D_manv_RT.gif';    	% File name for saved gif
iniName = 'SI_2D_manv_RT_ini.gif';	% File name for snapshot of initial condition 

figure
hold on
axis tight
axis equal
set(gcf,'Color',[1 1 1],'Position',[50 50 560*1.8 420*1.8])
set(gca,'XLim',[-1 12],'YLim',[-3 3],'FontSize',fs)
xlabel('x')
ylabel('y')

for ii=1:1:length(t)    % Go through time from 0 to tfinal
    h = zeros(2*n-3,1);     % handle for lines between each agent pair.
    kk = 1;
    jj = 1;
    plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'k.','MarkerSize',ms)    % Plot traj. dots
    jj = 6; % traj for leader
    plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'r.','MarkerSize',ms)    % Plot traj. dots
    % Line each connected agent pair
    for i0 = 1:n-1
        for j0 = i0+1:n
            if Adj(i0,j0) == 1
                h(kk) = line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                            [zz(ii,i0) zz(ii,j0)],'LineWidth',lw);
                kk = kk+1;
            end
        end
    end 
    % At time 0, draw initial configuration and creat initial gif
    if ii == 1
        for jj = 1:n
            hh1 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'rs','MarkerSize',ms);
        end
        legend(hh1,'Initial position','Location','Best')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im, map, gifName, 'GIF', 'WriteMode', 'overwrite',...
               'DelayTime', 0, 'LoopCount', inf);
        imwrite(im, map, iniName, 'GIF', 'WriteMode',...
                'overwrite', 'DelayTime', 0, 'LoopCount', inf);
%%%%% At time tfinal, draw final configuration%%%%%%%%%
    elseif ii == length(t) 
        for jj = 1:n
           hh2 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'ro','MarkerSize',ms);
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0) zz(ii,j0)],'LineWidth',lw);
                end
            end
        end 
        legend([hh1,hh2],'Initial position','final position',...
              'Location','NorthEast')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append',...
               'DelayTime', 5)
    else % in time between 0 and tfinal, make gif writemode be append
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append',...
               'DelayTime', 0)
    end
%     pause(pausetime) % Time pause between two trajectory do
    if (rem(ii,20)~=1)
        delete(h)
    end
end
h_msg = msgbox('Operation Completed', 'Success');
text(xx(length(t),:)'+0.15*ones(n,1),yy(length(t),:)'+0.1*ones(n,1),...
    zz(length(t),:)',['1';'2';'3';'4';'5';'6'],'FontSize',fs)
