%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates gif animation for single integrator model
% formation acquisition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_form_acq_results.mat'
lw = 2;                     % line width
ms = 12;                    % MakerSize
fs = 18;                    % Fontsize
pausetime = 0.1;            % Pause time
gifName = 'SI_acq.gif';     % File name for saved gif
iniName = 'SI_acq_ini.gif'; % File name for snapshot of initial condition 

figure
hold on
axis tight
set(gcf,'Color',[1 1 1],'Position',[50 95 560*1.8 420*1.8])
set(gca,'XLim',[-2 1],'YLim',[-1 1.6],'Box','on',...
    'DataAspectRatio',[1 1 1],'FontSize',fs,'YTick',[-1 -0.5 0 0.5 1 1.5])
xlabel('x')
ylabel('y')
for ii=1:2:length(t)    	% Go through time from 0 to tfinal
    h = zeros(2*n-3,1);     % handle for lines between each agent pair.
    kk = 1;
    plot(xx(ii,1:n),yy(ii,1:n),'k.','MarkerSize',ms)    % Plot traj. dots
    % Line each connected agent pair
    for i0 = 1:n-1
        for j0 = i0+1:n
            if Adj(i0,j0) == 1
                h(kk) = line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                            'LineWidth',lw);
                kk = kk+1;
            end
        end
    end 
    % At time 0, draw initial configuration and creat initial gif
    if ii == 1
        for jj = 1:n
            hh1 = plot(xx(ii,jj),yy(ii,jj),'rs','MarkerSize',ms);
        end
        legend(hh1,'Initial position','Location','NorthWest')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im, map, gifName, 'GIF', 'WriteMode', 'overwrite',...
               'DelayTime', 0, 'LoopCount', inf);
        imwrite(im, map, iniName, 'GIF', 'WriteMode',...
                'overwrite', 'DelayTime', 0, 'LoopCount', inf);
    elseif ii == length(t) % At time tfinal, draw final configuration
        for jj = 1:n
            hh2 = plot(xx(ii,jj),yy(ii,jj),'ro','MarkerSize',ms);
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        'LineWidth',lw);
                end
            end
        end
        legend([hh1,hh2],'initial position','final position',...
               'Location','NorthWest')
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