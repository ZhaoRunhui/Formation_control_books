%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 5 agents with single 
% integrator model performing formation acquisition task. Initial
% conditions are defined in SI_form_acq_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_form_acq_results.mat'
fs = 20;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 1;         % Linewidth 

%% Plot desired formation (convex pentagon)%%
figure
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
    'XTick',[-1 -0.5 0 0.5 1],'YTick',[-1 -0.5 0 0.5 1],...
    'XLim',[-1.2,1.2],'YLim',[-1.2,1.2])
plot(x_coor,y_coor,'-bo','MarkerSize',avgMs,'LineWidth',lw)
% Plot desired formation according to adjacency matrix
for ii = 1:n-1
    for jj = ii+1:n
        if Adj(ii,jj) == 1
            line([x_coor(ii) x_coor(jj)],[y_coor(ii),y_coor(jj)],...
                'LineWidth',lw)
        end
    end
end
grid on
xlabel('x')
ylabel('y')

%% Plot trajectory
figure
hold on
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs,...
    'YTick',[-1 -0.5 0 0.5 1])
for ii=1:1:length(t)
    for jj = 1:n
        % Plot trajectory dots
        plot(xx(ii,jj),yy(ii,jj),'k.','MarkerSize',avgMs) 
    end
    if ii == 1
        for jj = 1:n
            h1 = plot(xx(ii,jj),yy(ii,jj),'rs','MarkerSize',1.4*avgMs);
        end
    end
    if ii == length(t)
        for jj = 1:n
            h2 = plot(xx(ii,jj),yy(ii,jj),'bo','MarkerSize',1.4*avgMs);
            
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1
                    line([xx(ii,i0) xx(ii,j0)],[yy(ii,i0),yy(ii,j0)],...
                        'LineWidth',lw)
                end
            end
        end     
    end
    for jj = 1:n
        plot(xx(ii,jj),yy(ii,jj),'b.');
    end
end
legend([h1,h2],'Initial position','Final position','Location','NorthWest')
xlabel('x')
ylabel('y')
%% Plot distance error %%
nn = 1:10:length(t); % Plot distance error (e) using less points
figure
subplot(211)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e23(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e34(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e45(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e15(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
legend('e_1_2','e_2_3','e_3_4','e_4_5','e_1_5','Location','NorthEast')
xlabel('Time')
ylabel('e_i_j')
xlim([0 2.5])
subplot(212)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e13(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e14(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e24(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e25(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e35(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
legend('e_1_3','e_1_4','e_2_4','e_2_5','e_3_5','Location','NorthEast')
xlabel('Time')
ylabel('e_i_j')
xlim([0 2.5])

%% Plot control input %%
figure
subplot(211)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),u(1,nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),u(3,nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),u(5,nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),u(7,nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),u(9,nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
xlim([0 2.5])
xlabel('Time')
ylabel('u_i_x')
legend('u_1_x','u_2_x','u_3_x','u_4_x','u_5_x',...
       'Location','NorthEastOutside')
subplot(212)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),u(2,nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),u(4,nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),u(6,nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),u(8,nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),u(10,nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
xlim([0 2.5])
legend('u_1_y','u_2_y','u_3_y','u_4_y','u_5_y',...
       'Location','NorthEastOutside')
xlabel('Time')
ylabel('u_i_y')