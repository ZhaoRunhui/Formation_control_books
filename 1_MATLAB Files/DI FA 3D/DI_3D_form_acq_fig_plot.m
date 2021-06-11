%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 8 agents with double 
% integrator model performing formation acquisition task. Initial
% conditions is defined in DI_3D_form_acq_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'DI_3D_form_acq_results.mat'
fs = 14;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 1;         % Linewidth 

%% Plot desired formation (convex pentagon)%%
figure
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
   'XLim',[-2,2],'YLim',[-2,2],'ZLim',[-2,2])
for ii = 1:length(x_coor)
    plot3(x_coor(ii),y_coor(ii),z_coor(ii),'-ro','MarkerSize',avgMs,...
         'LineWidth',lw)
end
% Plot desired formation according to adjacency matrix
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
hold on
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
    'XLim',[-2,2],'YLim',[-2,2],'ZLim',[-2,2])
for ii=1:1:length(t)
    for jj = 1:n
        % Plot trajectory dots
        plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'b.','MarkerSize',avgMs) 
    end
    if ii == 1
        for jj = 1:n
            h1 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'rs',...
                       'MarkerSize',1.4*avgMs);
        end
    end
    if ii == length(t)
        for jj = 1:n
            h2 = plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'bo',...
                      'MarkerSize',1.4*avgMs);
            
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        [zz(ii,i0),zz(ii,j0)],'LineWidth',lw)
                end
            end
        end     
    end
    for jj = 1:n
        plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'b.');
    end
end
legend([h1,h2],'Initial position','Final position')
grid on 
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-2.5 2.5])
ylim([-2.5 2.5])
zlim([-2.5 2.5])
view([4,7,9])

%% Plot distance error %%
nn = 1:1:length(t); % Plot distance error (e) and velocity tracking error 
                    % (s) using less points
% The following code is only for 8 agents. Modification is needed for
% different number of agents
% figure one
figure
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
subplot(2,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e13(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e14(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e15(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e16(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
plot(t(nn),e17(nn),'-s','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 1]);
plot(t(nn),e18(nn),'-d','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[1 0.5 0]);
legend('e_1_2','e_1_3','e_1_4','e_1_5','e_1_6','e_1_7','e_1_8',...
     'Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% figure two
subplot(2,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e23(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e24(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e25(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e26(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e27(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
plot(t(nn),e28(nn),'-s','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 1]);
plot(t(nn),e34(nn),'-d','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[1 0.5 0]);
legend('e_2_3','e_2_4','e_2_5','e_2_6','e_2_7','e_2_8','e_3_4',...
     'Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% figure three
figure
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
subplot(2,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e35(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e36(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e37(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e38(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e45(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
plot(t(nn),e46(nn),'-s','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 1]);
plot(t(nn),e47(nn),'-d','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[1 0.5 0]);
legend('e_3_5','e_3_6','e_3_7','e_3_8','e_4_5','e_4_6','e_4_7',...
     'Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% figure four
subplot(2,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e48(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e56(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e57(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e58(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e67(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
plot(t(nn),e68(nn),'-s','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 1]);
plot(t(nn),e78(nn),'-d','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[1 0.5 0]);
legend('e_4_8','e_5_6','e_5_7','e_5_8','e_6_7','e_6_8','e_7_8',...
     'Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% Plot errors again
figure
subplot(211)
plot(t(nn),e12(nn),t(nn),e13(nn),t(nn),e14(nn),t(nn),e15(nn),...
     t(nn),e16(nn),t(nn),e17(nn),t(nn),e18(nn),t(nn),e23(nn),...
     t(nn),e24(nn),t(nn),e25(nn),t(nn),e26(nn),t(nn),e27(nn),...
     t(nn),e28(nn),t(nn),e34(nn),'LineWidth',1);
xlabel('Time')
ylabel('e_i_j')
grid

subplot(212)
plot(t(nn),e35(nn),t(nn),e36(nn),t(nn),e37(nn),t(nn),e38(nn),...
     t(nn),e45(nn),t(nn),e46(nn),t(nn),e47(nn),t(nn),e48(nn),...
     t(nn),e56(nn),t(nn),e57(nn),t(nn),e58(nn),t(nn),e67(nn),...
     t(nn),e68(nn),t(nn),e78(nn),'LineWidth',1);
xlabel('Time')
ylabel('e_i_j')
grid

%% Plot control input %%
figure
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
subplot(3,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),u( 1,nn),'LineWidth',1);
plot(t(nn),u( 4,nn),'LineWidth',1);
plot(t(nn),u( 7,nn),'LineWidth',1);
plot(t(nn),u(10,nn),'LineWidth',1);
plot(t(nn),u(13,nn),'LineWidth',1);
plot(t(nn),u(16,nn),'LineWidth',1);
plot(t(nn),u(19,nn),'LineWidth',1);
plot(t(nn),u(22,nn),'LineWidth',1);
xlabel('Time')
ylabel('u_i_x')

subplot(3,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),u( 2,nn),'LineWidth',1);
plot(t(nn),u( 5,nn),'LineWidth',1);
plot(t(nn),u( 8,nn),'LineWidth',1);
plot(t(nn),u(11,nn),'LineWidth',1);
plot(t(nn),u(14,nn),'LineWidth',1);
plot(t(nn),u(17,nn),'LineWidth',1);
plot(t(nn),u(20,nn),'LineWidth',1);
plot(t(nn),u(23,nn),'LineWidth',1);
xlabel('Time')
ylabel('u_i_y')

subplot(3,1,3)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),u( 3,nn),'LineWidth',1);
plot(t(nn),u( 6,nn),'LineWidth',1);
plot(t(nn),u( 9,nn),'LineWidth',1);
plot(t(nn),u(12,nn),'LineWidth',1);
plot(t(nn),u(15,nn),'LineWidth',1);
plot(t(nn),u(18,nn),'LineWidth',1);
plot(t(nn),u(21,nn),'LineWidth',1);
plot(t(nn),u(24,nn),'LineWidth',1);
xlabel('Time')
ylabel('u_i_z')

%% Plot s
s = (qv(:,(3*n+1):6*n))'-vf;

figure
set(gcf,'Color',[1 1 1],'Position',[150 50 560*1.8 480*1.8])
subplot(3,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),s( 1,nn),'LineWidth',1);
plot(t(nn),s( 4,nn),'LineWidth',1);
plot(t(nn),s( 7,nn),'LineWidth',1);
plot(t(nn),s(10,nn),'LineWidth',1);
plot(t(nn),s(13,nn),'LineWidth',1);
plot(t(nn),s(16,nn),'LineWidth',1);
plot(t(nn),s(19,nn),'LineWidth',1);
plot(t(nn),s(22,nn),'LineWidth',1);
xlabel('Time')
ylabel('s_i_x')

subplot(3,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),s( 2,nn),'LineWidth',1);
plot(t(nn),s( 5,nn),'LineWidth',1);
plot(t(nn),s( 8,nn),'LineWidth',1);
plot(t(nn),s(11,nn),'LineWidth',1);
plot(t(nn),s(14,nn),'LineWidth',1);
plot(t(nn),s(17,nn),'LineWidth',1);
plot(t(nn),s(20,nn),'LineWidth',1);
plot(t(nn),s(23,nn),'LineWidth',1);
xlabel('Time')
ylabel('s_i_y')

subplot(3,1,3)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),s( 3,nn),'LineWidth',1);
plot(t(nn),s( 6,nn),'LineWidth',1);
plot(t(nn),s( 9,nn),'LineWidth',1);
plot(t(nn),s(12,nn),'LineWidth',1);
plot(t(nn),s(15,nn),'LineWidth',1);
plot(t(nn),s(18,nn),'LineWidth',1);
plot(t(nn),s(21,nn),'LineWidth',1);
plot(t(nn),s(24,nn),'LineWidth',1);
xlabel('Time')
ylabel('s_i_z')