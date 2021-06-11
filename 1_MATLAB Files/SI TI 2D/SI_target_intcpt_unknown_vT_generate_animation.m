%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates gif animation for single integrator model
% target interception with unknown target velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_target_intcpt_unknown_vT_results.mat'
lw = 2;                                  % line width
ms = 12;                                 % MakerSize
fs = 24;                                 % fontsize
pausetime = 0.1;                         % Pause time
gifName = 'SI_intcpt_unknown_vT.gif';    % File name for saved gif
iniName = 'SI_intcpt_unknown_vT_ini.gif';% File name for snapshot of 
                                         % initial condition

figure
set(gcf,'Color',[1 1 1],'Position',[3 259 1440 408])
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',14,...
    'XLim',[-1 14],'YLim',[-2.1 2.1])
xlabel('x')
ylabel('y')
hold on
for ii=1:1:length(t)
    h = zeros(9,1);    
    kk = 1;
    % plot followers as dots
    plot(xx(ii,1:n-1),yy(ii,1:n-1),'k.','MarkerSize',ms)
    % plot target
    plot(xxt(ii),yyt(ii),'r+','MarkerSize',ms);
    % plot leader
    plot(xx(ii,n),yy(ii,n),'kh','MarkerSize',ms);
    for i0 = 1:n-1
        for j0 = i0+1:n
            if Adj(i0,j0) == 1
                h(kk) = line([xx(ii,i0) xx(ii,j0)],...
                             [yy(ii,i0),yy(ii,j0)],...
                            'LineWidth',lw);
                kk = kk+1;
            end
        end
    end 
    if ii == 1
        for jj = 1:n
            hh1 = plot(xx(ii,jj),yy(ii,jj),'rs','MarkerSize',ms);
        end
        hh3 = plot(xx(ii,n),yy(ii,n),'kh','MarkerSize',ms);
        hh4 = plot(xxt(ii),yyt(ii),'r+','MarkerSize',ms);
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        'LineWidth',lw);
                end
            end
        end 
        legend([hh1,hh3,hh4],'Followers initial pos.','Leader','Target')
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
        legend([hh1,hh2,hh3,hh4],'Follower initial pos.',...
                                 'Follower final pos.','Leader','Target')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append','DelayTime', 5)
    else	% in time between 0 and tfinal, make gif writemode be append
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append','DelayTime', 0)
    end
    pause(pausetime)
    delete(h)
end