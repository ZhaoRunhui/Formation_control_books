%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 6 agents with single 
% integrator model performing target interception with unknown target
% velocity. Initial conditions is defined in
% SI_target_intcpt_unknown_vT_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_target_intcpt_unknown_vT_results.mat'
fs = 24;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 2;         % Linewidth 

%% Plot desired formation (convex pentagon)%%
% figure
% hold on
% set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
%     'XTick',[-1 -0.5 0 0.5 1],'YTick',[-1 -0.5 0 0.5 1],...
%     'XLim',[-1.2,1.2],'YLim',[-1.2,1.2])
% plot(x_coor,y_coor,'bo','MarkerSize',avgMs)
% % Plot desired formation according to adjacency matrix
% for ii = 1:n-1
%     for jj = ii+1:n
%         if Adj(ii,jj) == 1
%             line([x_coor(ii) x_coor(jj)],[y_coor(ii),y_coor(jj)],...
%                 'LineWidth',lw)
%         else
%             line([x_coor(ii) x_coor(jj)],[y_coor(ii),y_coor(jj)],...
%                 'LineStyle','--','Color','r','LineWidth',lw)
%         end
%     end
% end
% grid on
% xlabel('x')
% ylabel('y')

%% Plot trajectory
figure
xlim([-1.2,12.5])
ylim([-2.7,1.5])
grid
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs)
for ii=1:1:length(t)
    hh1 = plot(xxt(ii),yyt(ii),'r+','MarkerSize',1.4*avgMs);
    hh2 = plot(xx(ii,n),yy(ii,n),'kh','MarkerSize',1.2*avgMs);
    if (ii==1)||(ii==16)||(ii==41)||(ii==71)||(ii==length(t))
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                         'LineWidth',lw)
                end
            end
        end
    end
%     pause(0.05) % Uncomment this if animation is needed
end
for i1 = 1:n-1
    hh3 = plot(xx(1,i1),yy(1,i1),'rs','MarkerSize',1.4*avgMs);
    hh4 = plot(xx(length(t),i1),yy(length(t),i1),'ro','MarkerSize',...
               1.4*avgMs);
end
legend([hh3,hh4,hh2,hh1],'Follower initial position',...
                         'Follower final position','Leader','Target')
xlabel('x')
ylabel('y')

%% Plot error and distance 
%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot error function 
%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = 1:1:length(t); % Plot error and distance using less points
figure
subplot(311)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'-+','LineWidth',1,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e23(nn),'-o','LineWidth',1,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e34(nn),'-*','LineWidth',1,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e45(nn),'.-','LineWidth',1,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e15(nn),'-x','LineWidth',1,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
legend('e_1_2','e_2_3','e_3_4','e_4_5','e_1_5','Location','NorthEast')
xlabel('Time')
ylabel('e_i_j')
xlim([0 9])
grid

subplot(312)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e13(nn),'-+','LineWidth',1,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e14(nn),'-o','LineWidth',1,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e24(nn),'-*','LineWidth',1,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e25(nn),'.-','LineWidth',1,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e35(nn),'-x','LineWidth',1,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
legend('e_1_3','e_1_4','e_2_4','e_2_5','e_3_5','Location','NorthEast')
xlabel('Time')
ylabel('e_i_j')
grid

subplot(313)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e16(nn),'-+','LineWidth',1,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e26(nn),'-o','LineWidth',1,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e36(nn),'-*','LineWidth',1,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e46(nn),'.-','LineWidth',1,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e56(nn),'-x','LineWidth',1,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
legend('e_1_6','e_2_6','e_3_6','e_4_6','e_5_6','Location','NorthEast')
xlabel('Time')
ylabel('e_i_j')
grid

%% Plot the control input
figure
subplot(211)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),u(1,nn),'LineWidth',1) 
plot(t(nn),u(3,nn),'LineWidth',1);
plot(t(nn),u(5,nn),'LineWidth',1)
plot(t(nn),u(7,nn),'LineWidth',1)
plot(t(nn),u(9,nn),'LineWidth',1)
plot(t(nn),u(11,nn),'k','LineWidth',1)
grid
xlabel('Time')
ylabel('u_i_x')

subplot(212)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),u(2,nn),'LineWidth',1) 
plot(t(nn),u(4,nn),'LineWidth',1) 
plot(t(nn),u(6,nn),'LineWidth',1) 
plot(t(nn),u(8,nn),'LineWidth',1) 
plot(t(nn),u(10,nn),'LineWidth',1) 
plot(t(nn),u(12,nn),'k','LineWidth',1)
grid
xlabel('Time')
ylabel('u_i_y')

% Plot errors again
figure
plot(t(nn),e12(nn),t(nn),e16(nn),t(nn),e23(nn),t(nn),e26(nn),...
     t(nn),e34(nn),t(nn),e36(nn),t(nn),e45(nn),t(nn),e46(nn),...
     t(nn),e56(nn),t(nn),e15(nn),t(nn),e13(nn),t(nn),e14(nn),...
     t(nn),e24(nn),t(nn),e25(nn),t(nn),e35(nn),'LineWidth',1);
xlabel('Time')
ylabel('e_i_j')
grid