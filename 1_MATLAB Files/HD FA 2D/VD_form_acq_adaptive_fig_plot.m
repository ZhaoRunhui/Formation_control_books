%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 5 agents with vehicle 
% dynamics model performing formation acquisition task. Initial conditions
% are defined in VD_form_acq_adaptive_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
% Subfunction: plot_vehicle.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'VD_form_acq_adaptive_results.mat'
fs = 24;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 2;         % Linewidth 

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
xlim([-1.05,1])
ylim([-0.9,1.35])
set(gca,'DataAspectRatio',[1 1 1],'Box','on','FontSize',fs)
for ii=1:1:length(t)
    for jj = 1:n
        plot(qhx(ii,jj),qhy(ii,jj),'k.','MarkerSize',avgMs)
        if ii == 1
           h1 = plot(qhx(ii,jj),qhy(ii,jj),'rs','MarkerSize',1.4*avgMs,...
                     'LineWidth',lw);
        end
        if ii == length(t)
           h2 = plot(qhx(ii,jj),qhy(ii,jj),'ro','MarkerSize',1.4*avgMs,...
                     'LineWidth',lw);
        end
    end
    if (ii==1)||(ii==length(t))
        for jj = 1:n
            plot_vehicle(qcx(ii,jj),qcy(ii,jj),theta(ii,jj),L(jj));
        end
        for i0 = 1:n-1
            for j0 = i0+1:n
                if Adj(i0,j0) == 1 && (ii==length(t))
                   line([qhx(ii,i0) qhx(ii,j0)],[qhy(ii,i0),qhy(ii,j0)],...
                    'LineWidth',lw)
                end
            end
        end
    end
%     pause(0.05) % Uncomment this if animation is needed
end
legend([h1,h2],'Initial position','Final position','Location','NorthEast')
xlabel('x (m)')
ylabel('y (m)')
%% Plot distance error %%
nn = 1:8:length(t); % Plot distance error (e) and velocity tracking error 
                    % (s) using less points
figure
subplot(211)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e23(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e34(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e45(nn),'.-','LineWidth',lw,'MarkerSize',2.0*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e15(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
xlabel('Time (s)')
ylabel('${e}_{ij}\mbox{ (m)}$','Interpreter','latex')
h1 = legend('${e}_{12}$','${e}_{23}$','${e}_{34}$',...
        '${e}_{45}$','${e}_{15}$','Location','SouthEast');
set(h1,'Interpreter','latex')
xlim([0 10])
ylim([-0.4 0.1])
subplot(212)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e13(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[0 0 1]);
plot(t(nn),e14(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,...
     'Color',[0 0.5 0]);
plot(t(nn),e24(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,...
     'Color',[1 0 0]);
plot(t(nn),e25(nn),'.-','LineWidth',lw,'MarkerSize',2.0*avgMs,...
     'Color',[0 0 0]);
plot(t(nn),e35(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,...
     'Color',[0.75 0 0.75]);
xlabel('Time (s)')
ylabel('${e}_{ij}\mbox{ (m)}$','Interpreter','latex')
h2 = legend('${e}_{13}$','${e}_{14}$','${e}_{24}$',...
        '${e}_{25}$','${e}_{35}$','Location','SouthEast');
set(h2,'Interpreter','latex')
xlim([0 10])
ylim([-0.6 0.15])


%% Plot the actual control input %%
figure
subplot(211)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),ubar(1,nn),'LineWidth',1);
plot(t(nn),ubar(3,nn),'LineWidth',1);
plot(t(nn),ubar(5,nn),'LineWidth',1);
plot(t(nn),ubar(7,nn),'LineWidth',1);
plot(t(nn),ubar(9,nn),'LineWidth',1);
xlim([0 10])
ylim([-2.5 2])
xlabel('Time (s)')
ylabel('$\bar{u}_{ix}\mbox{ (N)}$','Interpreter','latex')
grid

subplot(212)
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),ubar(2,nn),'LineWidth',1);
plot(t(nn),ubar(4,nn),'LineWidth',1);
plot(t(nn),ubar(6,nn),'LineWidth',1);
plot(t(nn),ubar(8,nn),'LineWidth',1);
plot(t(nn),ubar(10,nn),'LineWidth',1);
xlim([0 10])
ylim([-0.35 0.35])
xlabel('Time (s)')
ylabel('$\bar{u}_{iy}\mbox{ (N-m)}$','Interpreter','latex')
grid

%% Plot error of estimation dynamics
for ii = 1:n
    figure(ii+4)
    hold on
    set(gca,'Box','on','FontSize',fs)
    plot(t(nn),phihat(nn,6*ii-5),'-+','LineWidth',lw,...
        'MarkerSize',1.2*avgMs,'Color',[0 0 1]);
    plot(t(nn),phihat(nn,6*ii-4),'-o','LineWidth',lw,...
        'MarkerSize',0.8*avgMs,'Color',[0 0.5 0]);
    plot(t(nn),phihat(nn,6*ii-3),'-*','LineWidth',lw,...
        'MarkerSize',1.2*avgMs,'Color',[1 0 0]);
    plot(t(nn),phihat(nn,6*ii-2),'.-','LineWidth',lw,...
        'MarkerSize',2*avgMs,'Color',[0 0 0]);
    plot(t(nn),phihat(nn,6*ii-1),'-x','LineWidth',lw,...
        'MarkerSize',1.4*avgMs,'Color',[0.75 0 0.75]);
    plot(t(nn),phihat(nn,6*ii),'--','LineWidth',lw,'Color',[0 0 0]);
    xlabel('Time (s)')
    ylabel(['$\hat{\phi}_{',num2str(ii),'}$'],'Interpreter','latex')
    h = legend(['$[\hat{\phi}_{',num2str(ii),'}]_{1}$'],...
               ['$[\hat{\phi}_{',num2str(ii),'}]_{2}$'],...
               ['$[\hat{\phi}_{',num2str(ii),'}]_{3}$'],...
               ['$[\hat{\phi}_{',num2str(ii),'}]_{4}$'],...
               ['$[\hat{\phi}_{',num2str(ii),'}]_{5}$'],...
               ['$[\hat{\phi}_{',num2str(ii),'}]_{6}$'],...
               'Location','NorthEast');
    set(h,'Interpreter','latex')
end

% Plot errors again
figure
plot(t(nn),e12(nn),t(nn),e23(nn),t(nn),e34(nn),t(nn),e45(nn),...
     t(nn),e15(nn),t(nn),e13(nn),t(nn),e14(nn),t(nn),e24(nn),...
     t(nn),e25(nn),t(nn),e35(nn),'LineWidth',1);
xlabel('Time')
ylabel('e_i_j')
grid