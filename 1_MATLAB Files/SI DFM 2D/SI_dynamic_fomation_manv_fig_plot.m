%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 5 agents with single 
% integrator model performing dynamic formation maneuvering task. Initial
% conditions are defined in SI_dynamic_fomation_manv_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_dynamic_fomation_manv_results.mat'
fs = 24;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 2;         % Linewidth 
tf = max(t);

%% Plot trajectory
figure
xlim([-1,17]); %22])
ylim([-2.1,2])
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs)
for ii=1:1:length(t)
    for jj = 1:n
        plot(xx(ii,jj),yy(ii,jj),'k.','MarkerSize',avgMs);
    end
    if (ii==1)||(ii==41)||(ii==81)||(ii==101)||...
            (ii==121)||(ii==length(t))
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
for i1 = 1:n
    hh1 = plot(xx(1,i1),yy(1,i1),'rs','MarkerSize',1.4*avgMs);
    hh2 = plot(xx(length(t),i1),yy(length(t),i1),'ro',...
               'MarkerSize',1.4*avgMs);
end
rectangle('Position',[10.4,0.5,2.2,0.3],'Facecolor','blue',...
          'Edgecolor','none')
rectangle('Position',[10.4,-1.3,2.2,0.3],'Facecolor','blue',...
          'Edgecolor','none')
legend([hh1,hh2],'Intitial position', 'Final position')
xlabel('x')
ylabel('y')

%% Plot distance error %%
nn = 1:1:length(t); % Plot distance error (e) using less points
figure
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'LineWidth',lw,'Color',[0 0 1]);
plot(t(nn),e23(nn),'LineWidth',lw,'Color',[0 0.5 0]);
plot(t(nn),e34(nn),'LineWidth',lw,'Color',[1 0 0]);
plot(t(nn),e45(nn),'LineWidth',lw,'Color',[0 0 0]);
plot(t(nn),e15(nn),'LineWidth',lw,'Color',[0.75 0 0.75]);
plot(t(nn),e13(nn),'LineWidth',lw,'Color',[0 0 1]);
plot(t(nn),e14(nn),'LineWidth',lw,'Color',[0 0.5 0]);
plot(t(nn),e24(nn),'LineWidth',lw,'Color',[1 0 0]);
plot(t(nn),e25(nn),'LineWidth',lw,'Color',[0 0 0]);
plot(t(nn),e35(nn),'LineWidth',lw,'Color',[0.75 0 0.75]);
xlabel('Time')
ylabel('e_i_j')
xlim([0 1])



%% Plot control input %%
figure
subplot(211)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),u(1,nn),'LineWidth',lw,'Color',[0 0 1]);
plot(t(nn),u(3,nn),'LineWidth',lw,'Color',[0 0.5 0]);
plot(t(nn),u(5,nn),'LineWidth',lw,'Color',[1 0 0]);
plot(t(nn),u(7,nn),'LineWidth',lw,'Color',[0 0 0]);
plot(t(nn),u(9,nn),'LineWidth',lw,'Color',[0.75 0 0.75]);
xlim([0 tf])
xlabel('Time')
ylabel('u_i_x')
legend('u_1_x','u_2_x','u_3_x','u_4_x','u_5_x','Location','NorthEastOutside')
subplot(212)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),u(2,nn),'LineWidth',lw,'Color',[0 0 1]);
plot(t(nn),u(4,nn),'LineWidth',lw,'Color',[0 0.5 0]);
plot(t(nn),u(6,nn),'LineWidth',lw,'Color',[1 0 0]);
plot(t(nn),u(8,nn),'LineWidth',lw,'Color',[0 0 0]);
plot(t(nn),u(10,nn),'LineWidth',lw,'Color',[0.75 0 0.75]);
xlim([0 tf])
legend('u_1_y','u_2_y','u_3_y','u_4_y','u_5_y','Location','NorthEastOutside')
xlabel('Time')
ylabel('u_i_y')

