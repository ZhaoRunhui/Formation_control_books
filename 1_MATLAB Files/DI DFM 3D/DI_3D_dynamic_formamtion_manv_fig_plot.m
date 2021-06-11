%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plot the result from
%   DI_3D_dynamic_formamtion_manv_main.m
%   
% output:
%   1. desired formation
%   2. trajectory of the righd body
%   3. distance error of edges
%   4. control inputs
%   5. velocity tracking error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'DI_3D_dynamic_formamtion_manv_results.mat'
fs = 14;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 2;         % Linewidth 

%% Plot desired formation %%
figure
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
hold on

set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
    'XLim',[-2,2],'YLim',[-2,2],'ZLim',[-2,2])
for ii = 1:length(x_coor)
    plot3(x_coor(ii),y_coor(ii),z_coor(ii),'-ro','MarkerSize',avgMs,...
          'LineWidth',lw)
    % Plot desired formation according to adjacency matrix
end
for ii = 1:n-1
    for jj = ii+1:n
        if Adj(ii,jj) == 1
            line([x_coor(ii) x_coor(jj)],[y_coor(ii),y_coor(jj)],...
                 [z_coor(ii),z_coor(jj)],'LineWidth',lw);
         
        end
    end
end
grid on
view([4,7,9])
xlabel('x')
ylabel('y')
zlabel('z')

%% Plot trajectory
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs)
grid on
for ii=1:2:length(t)
   %%% plot the trajectory of agent 9
   jj = 9;
   plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'k.','MarkerSize',avgMs);
   %%% snapshot setup %%
    if (ii==1)||(ii==41)||(ii==81)||(ii==121)||(ii==161)||(ii== 201)...
            ||(ii==241)||(ii==281)||(ii==length(t))
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0) zz(ii,j0)],'LineWidth',1)
                end
            end
        end
    end
%     pause(0.05) % Uncomment this if animation is needed
end
for i1 = 1:n
    hh1 = plot3(xx(1,i1),yy(1,i1),zz(1,i1),'rs','MarkerSize',1.4*avgMs);
    hh2 = plot3(xx(length(t),i1),yy(length(t),i1),zz(length(t),i1),'ro',...
               'MarkerSize',1.4*avgMs);
end
legend([hh1,hh2],'Intitial position', 'Final position')
xlabel('x')
ylabel('y')
zlabel('z')
view([10,-3,6])

%% Plot distance error %%
nn = 1:1:length(t); % Plot distance error (e), control input (u) and  
                    % velocity tracking error (s)
% figure one
lw=1;
figure
subplot(2,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'LineWidth',lw);
plot(t(nn),e13(nn),'LineWidth',lw);
plot(t(nn),e14(nn),'LineWidth',lw);
plot(t(nn),e15(nn),'LineWidth',lw);
plot(t(nn),e16(nn),'LineWidth',lw);
plot(t(nn),e17(nn),'LineWidth',lw);
plot(t(nn),e18(nn),'LineWidth',lw);
plot(t(nn),e19(nn),'LineWidth',lw);
plot(t(nn),e23(nn),'LineWidth',lw);
plot(t(nn),e24(nn),'LineWidth',lw);
plot(t(nn),e25(nn),'LineWidth',lw);
plot(t(nn),e26(nn),'LineWidth',lw);
plot(t(nn),e27(nn),'LineWidth',lw);
plot(t(nn),e28(nn),'LineWidth',lw);
plot(t(nn),e29(nn),'LineWidth',lw);
plot(t(nn),e34(nn),'LineWidth',lw);
plot(t(nn),e35(nn),'LineWidth',lw);
plot(t(nn),e36(nn),'LineWidth',lw);
xlabel('Time')
ylabel('e_i_j')
xlim([0 6])
ylim([-0.8,1.2])

% figure two
subplot(2,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e37(nn),'LineWidth',lw);
plot(t(nn),e38(nn),'LineWidth',lw);
plot(t(nn),e39(nn),'LineWidth',lw);
plot(t(nn),e45(nn),'LineWidth',lw);
plot(t(nn),e46(nn),'LineWidth',lw);
plot(t(nn),e47(nn),'LineWidth',lw);
plot(t(nn),e48(nn),'LineWidth',lw);
plot(t(nn),e49(nn),'LineWidth',lw);
plot(t(nn),e56(nn),'LineWidth',lw);
plot(t(nn),e57(nn),'LineWidth',lw);
plot(t(nn),e58(nn),'LineWidth',lw);
plot(t(nn),e59(nn),'LineWidth',lw);
plot(t(nn),e67(nn),'LineWidth',lw);
plot(t(nn),e68(nn),'LineWidth',lw);
plot(t(nn),e69(nn),'LineWidth',lw);
plot(t(nn),e78(nn),'LineWidth',lw);
plot(t(nn),e79(nn),'LineWidth',lw);
plot(t(nn),e89(nn),'LineWidth',lw);
xlabel('Time')
ylabel('e_i_j')
xlim([0 6])
ylim([-1.1,0.8])


%% Plot control input %%
figure
subplot(3,1,1)
grid on
hold on
plot(t(nn),u( 1,nn),'LineWidth',lw);
plot(t(nn),u( 4,nn),'LineWidth',lw);
plot(t(nn),u( 7,nn),'LineWidth',lw);
plot(t(nn),u(10,nn),'LineWidth',lw);
plot(t(nn),u(13,nn),'LineWidth',lw);
plot(t(nn),u(16,nn),'LineWidth',lw);
plot(t(nn),u(19,nn),'LineWidth',lw);
plot(t(nn),u(22,nn),'LineWidth',lw);
plot(t(nn),u(25,nn),'LineWidth',lw);
xlim([0 6])
xlabel('Time')
ylabel('u_i_x')

subplot(3,1,2)
grid on
hold on
plot(t(nn),u( 2,nn),'LineWidth',lw);
plot(t(nn),u( 5,nn),'LineWidth',lw);
plot(t(nn),u( 8,nn),'LineWidth',lw);
plot(t(nn),u(11,nn),'LineWidth',lw);
plot(t(nn),u(14,nn),'LineWidth',lw);
plot(t(nn),u(17,nn),'LineWidth',lw);
plot(t(nn),u(20,nn),'LineWidth',lw);
plot(t(nn),u(23,nn),'LineWidth',lw);
plot(t(nn),u(26,nn),'LineWidth',lw);
xlabel('Time')
ylabel('u_i_y')
xlim([0 6])

subplot(3,1,3)
grid on
hold on
plot(t(nn),u( 3,nn),'LineWidth',lw);
plot(t(nn),u( 6,nn),'LineWidth',lw);
plot(t(nn),u( 9,nn),'LineWidth',lw);
plot(t(nn),u(12,nn),'LineWidth',lw);
plot(t(nn),u(15,nn),'LineWidth',lw);
plot(t(nn),u(18,nn),'LineWidth',lw);
plot(t(nn),u(21,nn),'LineWidth',lw);
plot(t(nn),u(24,nn),'LineWidth',lw);
plot(t(nn),u(27,nn),'LineWidth',lw);
xlabel('Time')
ylabel('u_i_z')
xlim([0 6])

%% Plot s
figure
subplot(3,1,1)
grid on
hold on
plot(t(nn),s( 1,nn),'LineWidth',lw);
plot(t(nn),s( 4,nn),'LineWidth',lw);
plot(t(nn),s( 7,nn),'LineWidth',lw);
plot(t(nn),s(10,nn),'LineWidth',lw);
plot(t(nn),s(13,nn),'LineWidth',lw);
plot(t(nn),s(16,nn),'LineWidth',lw);
plot(t(nn),s(19,nn),'LineWidth',lw);
plot(t(nn),s(22,nn),'LineWidth',lw);
plot(t(nn),s(25,nn),'LineWidth',lw);
xlim([0 6])
xlabel('Time')
ylabel('s_i_x')

subplot(3,1,2)
grid on
hold on
plot(t(nn),s( 2,nn),'LineWidth',lw);
plot(t(nn),s( 5,nn),'LineWidth',lw);
plot(t(nn),s( 8,nn),'LineWidth',lw);
plot(t(nn),s(11,nn),'LineWidth',lw);
plot(t(nn),s(14,nn),'LineWidth',lw);
plot(t(nn),s(17,nn),'LineWidth',lw);
plot(t(nn),s(20,nn),'LineWidth',lw);
plot(t(nn),s(23,nn),'LineWidth',lw);
plot(t(nn),s(26,nn),'LineWidth',lw);
xlabel('Time')
ylabel('s_i_y')
xlim([0 6])

subplot(3,1,3)
grid on
hold on
plot(t(nn),s( 3,nn),'LineWidth',lw);
plot(t(nn),s( 6,nn),'LineWidth',lw);
plot(t(nn),s( 9,nn),'LineWidth',lw);
plot(t(nn),s(12,nn),'LineWidth',lw);
plot(t(nn),s(15,nn),'LineWidth',lw);
plot(t(nn),s(18,nn),'LineWidth',lw);
plot(t(nn),s(21,nn),'LineWidth',lw);
plot(t(nn),s(24,nn),'LineWidth',lw);
plot(t(nn),s(27,nn),'LineWidth',lw);
xlabel('Time')
ylabel('s_i_z')
xlim([0 6])