%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates gif animation for 3D Double integrator model
% target interception with unknown target Acceleration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'DI_3D_target_intcpt_unknown_aT_results.mat'
lw = 2;                                     % line width
ms = 12;                                    % MakerSize
fs = 24;                                    % fontsize
pausetime = 0.1;                            % Pause time
gifName = 'DI_3D_intcpt_unknown_aT.gif';    % File name for saved gif
iniName = 'DI_3D_intcpt_unknown_aT_ini.gif';% File name for snapshot 
                                            % of initial condition
figure
hold on
axis tight
set(gcf,'Color',[1 1 1],'Position',[50 50 560*2.5 420*2])
set(gca,'XLim',[-3 3],'YLim',[-2 34],'ZLim',[-2.1,4],'Box','on',...
    'DataAspectRatio',[1 1 1],'FontSize',fs)
xlabel('x')
ylabel('y')
zlabel('z')
view([20,-9,8])

for ii=1:1:length(t)
    h = zeros(3*n-6,1);    
    kk = 1;
    % plot followers as dots
    plot3(xx(ii,1:n-1),yy(ii,1:n-1),zz(ii,1:n-1),'k.','MarkerSize',ms)
    % plot target
    plot3(xxt(ii),yyt(ii),zzt(ii),'r+','MarkerSize',ms);
    % plot leader
    plot3(xx(ii,n),yy(ii,n),zz(ii,n),'kh','MarkerSize',ms);
    for i0 = 1:n-1
        for j0 = i0+1:n
            if Adj(i0,j0) == 1
                h(kk) = line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                            [zz(ii,i0) zz(ii,j0)],'LineWidth',lw);
                kk = kk+1;
            end
        end
    end
    % at time 0, draw initial configuration and creat initial gif
    if ii == 1
        for jj = 1:n-1
            hh1 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),...
                       'rs','MarkerSize',ms); % followers
        end
        %leader
        hh3 = plot3(xx(ii,n),yy(ii,n),zz(ii,n),'kh','MarkerSize',ms);
        % target
        hh4 = plot3(xxt(ii),yyt(ii),zzt(ii),'r+','MarkerSize',ms); 
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0),zz(ii,j0)],'LineWidth',lw);
                end
            end
        end 
        legend([hh1,hh3,hh4],'Followers initial pos.','Leader','Target',...
            'Location','Best')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im, map, gifName, 'GIF', 'WriteMode', 'overwrite',...
                'DelayTime', 0, 'LoopCount', inf);
        imwrite(im, map, iniName, 'GIF', 'WriteMode',...
                'overwrite', 'DelayTime', 0, 'LoopCount', inf);
   % At time final, draw final configuration
    elseif ii == length(t)
        for jj = 1:n
            hh2 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'ro',...
                        'MarkerSize',ms);
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0),zz(ii,j0)],'LineWidth',lw);
                end
            end
        end 
        legend([hh1,hh2,hh3,hh4],'Follower initial pos.',...
                                 'Follower final pos.','Leader',...
                                  'Target','Location','Best')
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append','DelayTime', 5)
    else
        % in the time between 0 and tfinal, make gif writemode to be append
        f = getframe(gcf);
        f = frame2im(f);
        [im,map] = rgb2ind(f,128);
        imwrite(im,map,gifName,'gif','WriteMode','append','DelayTime', 0)
    end
%     pause(pausetime) % Time pause between two trajectory dots
    delete(h)
end
text(xx(length(t),:)'-0.05*ones(n,1),yy(length(t),:)'+0.05*ones(n,1),...
    zz(length(t),:)',['1';'2';'3';'4';'5';'6';'7';'8';'9'],...
    'FontSize',0.6*fs)
