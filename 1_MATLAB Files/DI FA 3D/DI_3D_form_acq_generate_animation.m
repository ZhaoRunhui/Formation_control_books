%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates gif animation for Double integrator model
% formation acquisition in 3D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'DI_3D_form_acq_results.mat'
lw = 2;                       % line width
ms = 12;                      % MakerSize
fs = 18;                      % Fontsize
pausetime = 0.1;              % Pause time
gifName = 'DI_3D_acq.gif';    % File name for saved gif
iniName = 'DI_3D_acq_ini.gif';% File name for snapshot of initial condition 

figure
hold on
axis tight
set(gcf,'Color',[1 1 1],'Position',[50 50 560*1.8 480*1.8])
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
    'XLim',[-2,2],'YLim',[-2,2],'ZLim',[-2,2])
xlabel('x')
ylabel('y')
zlabel('z')
view([14,7,10])

for ii=1:1:length(t)        % Iterate time from 0 to tfinal
    h = zeros(3*n-6,1);     % handle for lines between each agent pair.
    kk = 1;
for jj = 1:n;
    % Plot traj. dots
    plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'k.','MarkerSize',ms)  
end
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
    elseif ii == length(t) % At time tfinal, draw final configuration
        for jj = 1:n
           hh2 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'ro','MarkerSize',ms);
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0) zz(ii,j0)],'LineWidth',lw);
                    kk = kk+1;
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
    else % in time between 0 and tfinal, make gif writemode to be append
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append',...
               'DelayTime', 0)
    end
%     pause(pausetime) % Time pause between two trajectory dots
   delete(h)
end
text(xx(length(t),:)'-0.05*ones(n,1),yy(length(t),:)'+0.05*ones(n,1),...
    zz(length(t),:)',['1';'2';'3';'4';'5';'6';'7';'8'],'FontSize',fs)
