%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 9 agents of 3D Double 
% integrator model performing target interception task. Initial conditions
% are defined in DI_3D_target_intcpt_unknown_aT_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'DI_3D_target_intcpt_unknown_aT_results.mat'
fs = 14;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 2;         % Linewidth 

%% Plot desired formation (unit cube)%%
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
    'XTick',[0 0.5 1],'YTick',[0 0.5 1],...
    'ZTick',[0 0.5 1],'XLim',[-1.2,1.2],'YLim',[-1.2,1.2],...
    'ZLim',[-1.2,1.2])
 plot3(x_coor,y_coor,z_coor,'ro','MarkerSize',avgMs,'LineWidth',lw)
% Plot desired formation according to adjacency matrix
for ii = 1:n-1
    for jj = ii+1:n
        if Adj(ii,jj) == 1
            line([x_coor(ii),x_coor(jj)],[y_coor(ii),y_coor(jj)],...
                [z_coor(ii),z_coor(jj)],'LineWidth',lw)
        end
    end
end
grid on
xlabel('x')
ylabel('y')
zlabel('z')
text(x_coor-0.05*ones(n,1),y_coor+0.05*ones(n,1),...
    z_coor,['1';'2';'3';'4';'5';'6';'7';'8';'9'],'FontSize',fs)
view([10,7,4])

%% Plot trajectory
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs)
grid on
for ii=1:1:length(t)
    for jj = 1:n
    end
    hh1 = plot3(xxt(ii),yyt(ii),zzt(ii),'r+','MarkerSize',1.4*avgMs);
    hh2 = plot3(xx(ii,n),yy(ii,n),zz(ii,n),'kh','MarkerSize',1.2*avgMs);
    if (ii==1)||(ii==length(t))
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0) zz(ii,j0)],'LineWidth',lw)
                end
            end
        end
    end
%     pause(0.05) % Uncomment this if animation is needed
end
for i1 = 1:n
    hh3 = plot3(xx(1,i1),yy(1,i1),zz(1,i1),'rs','MarkerSize',1.4*avgMs);
    hh4 = plot3(xx(length(t),i1),yy(length(t),i1),...
                zz(length(t),i1),'ro','MarkerSize',1.4*avgMs);
end
legend([hh1,hh2,hh3,hh4],'Target', 'Leader',...
       'Follower intitial pos.','Follower final pos.')
xlabel('x')
ylabel('y')
zlabel('z')
text(xx(length(t),:)'-0.05*ones(n,1),yy(length(t),:)'+0.05*ones(n,1),...
     zz(length(t),:)',['1';'2';'3';'4';'5';'6';'7';'8';'9'],'FontSize',fs)
view([10,-3,6])

%% Plot distance error %%
nn = 1:2:length(t); % Plot distance error (e) using less points

% fig one
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
subplot(3,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e14(nn))
plot(t(nn),e15(nn))
plot(t(nn),e19(nn))
plot(t(nn),e23(nn))
plot(t(nn),e26(nn))
plot(t(nn),e29(nn))
plot(t(nn),e34(nn))
legend('e_1_4','e_1_5','e_1_9','e_2_3','e_2_6','e_2_9','e_3_4','Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% fig two
subplot(3,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e37(nn))
plot(t(nn),e39(nn))
plot(t(nn),e45(nn))
plot(t(nn),e48(nn))
plot(t(nn),e49(nn))
plot(t(nn),e56(nn))
plot(t(nn),e58(nn))
legend('e_3_7','e_3_9','e_4_5','e_4_8','e_4_9','e_5_6','e_5_8','Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% fig three
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
subplot(3,1,3)

grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e59(nn))
plot(t(nn),e67(nn))
plot(t(nn),e68(nn))
plot(t(nn),e69(nn))
plot(t(nn),e78(nn))
plot(t(nn),e79(nn))
plot(t(nn),e89(nn))
legend('e_5_9','e_6_7','e_6_8','e_6_9','e_7_8','e_7_9','e_8_9','Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')


% additional figure for distance between leader and followers
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])

grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e19(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e29(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e39(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e49(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e59(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
plot(t(nn),e69(nn),'-s','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 1]);
plot(t(nn),e79(nn),'-d','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[1 0.5 0]);
plot(t(nn),e89(nn),'-k','LineWidth',lw,'MarkerSize',0.6*avgMs,...
     'Color',[0 0.5 0.5]);
legend('e_1_9','e_2_9','e_3_9','e_4_9','e_5_9','e_6_9','e_7_9','e_8_9',...
      'Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% Plot errors again
lw=1;
figure
subplot(2,2,1)
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
xlim([0 10])

subplot(2,2,3)
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
xlim([0 10])
ylim([-0.7,0.5])

subplot(2,2,2)
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
xlim([10 30])

subplot(2,2,4)
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
xlim([10 30])

%% Plot control input %%
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
subplot(3,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),u( 1,nn),'LineWidth',lw);
plot(t(nn),u( 4,nn),'LineWidth',lw);
plot(t(nn),u( 7,nn),'LineWidth',lw);
plot(t(nn),u(10,nn),'LineWidth',lw);
plot(t(nn),u(13,nn),'LineWidth',lw);
plot(t(nn),u(16,nn),'LineWidth',lw);
plot(t(nn),u(19,nn),'LineWidth',lw);
plot(t(nn),u(22,nn),'LineWidth',lw);
plot(t(nn),u(25,nn),'LineWidth',lw);
xlim([0 20])
xlabel('Time')
ylabel('u_i_x')

subplot(3,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
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
xlim([0 20])

subplot(3,1,3)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
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
xlim([0 20])
h_msg = msgbox('Operation Completed', 'Success');

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
xlim([0 10])
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
xlim([0 10])

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
xlim([0 10])

