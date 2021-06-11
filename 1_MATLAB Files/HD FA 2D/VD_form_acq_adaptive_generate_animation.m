%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates gif animation for vehicle dynamic model
% formation acquisition with unknown dynamics.
% subfunction: plot_vehicle.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'VD_form_acq_adaptive_results.mat'
lw = 2;                              % line width
ms = 12;                             % MakerSize
fs = 24;                             % Fontsize
pausetime = 0.01;                    % Pause time
gifName = 'VD_form_acq_adap.gif';    % File name for saved gif
iniName = 'VD_form_acq_adap_ini.gif';% File name for snapshot of ini. cond.

figure
hold on
xlim([-1.05,1])
ylim([-0.9,1.35])
set(gcf,'Color',[1 1 1],'Position',[4 32 1671 950])
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs)
xlabel('x(m)')
ylabel('y(m)')
for ii=1:10:length(t)    % Go through time from 0 to tfinal
    h = zeros(7,1);
    hh = zeros(4,n);
    kk = 1;
    plot(qhx(ii,1:n),qhy(ii,1:n),'k.','MarkerSize',ms)
    for jj = 1:n
        hh(:,jj) = plot_vehicle(qcx(ii,jj),qcy(ii,jj),theta(ii,jj),L(jj));
    end
    for i0 = 1:n-1
        for j0 = i0+1:n
            if Adj(i0,j0) == 1
                h(kk) = line([qhx(ii,i0) qhx(ii,j0)],...
                            [qhy(ii,i0),qhy(ii,j0)],'LineWidth',lw);
                kk = kk+1;
            end
        end
    end 
    % At time 0, draw initial configuration and creat initial gif
    if ii == 1
        h1 = plot(qhx(ii,1:n),qhy(ii,1:n),'bs','MarkerSize',ms+2);
        legend(h1,'initial position','Location','NorthEast')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im, map, gifName, 'GIF', 'WriteMode', 'overwrite',...
               'DelayTime', 0, 'LoopCount', inf);
        imwrite(im, map, iniName, 'GIF', 'WriteMode',...
                'overwrite', 'DelayTime', 0, 'LoopCount', inf);
    elseif ii == length(t) % At time tfinal, draw final configuration
        h2 = plot(qhx(ii,1:n),qhy(ii,1:n),'bo','MarkerSize',ms+2);
        for jj = 1:n
            plot_vehicle(qcx(ii,jj),qcy(ii,jj),theta(ii,jj),L(jj));
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([qhx(ii,i0) qhx(ii,j0)],...
                         [qhy(ii,i0),qhy(ii,j0)],'LineWidth',lw);
                end
            end
        end
        legend([h1,h2],'initial position','final position',...
                        'Location','NorthEast')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append',...
               'DelayTime', 5)
    else 
        % in time between 0 and tfinal, make gif writemode to be append
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append',...
               'DelayTime', 0)
    end
    pause(pausetime) % Time pause between two trajectory dots
    delete(h,hh)
end