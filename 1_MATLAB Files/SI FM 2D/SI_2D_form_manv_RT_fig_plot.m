%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program plots publishable quality figures for 5 agents plus 1 leader
% with single integrator model performing formation maneuvering task. 
% Initial conditions is defined in SI_form_manv_RT_main.m
% Modification is required if any of the following is changed:
% initial condition, desired formation, number of agents, connecting of
% graph etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
load 'SI_2D_form_manv_RT_results.mat'
fs = 14;        % Font size in the figure
avgMs = 10;     % Average marker size
lw = 2;         % Linewidth 
n = 6;

%% Plot desired formation (unit cube)%%
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
hold on
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
axis equal
text(x_coor-0.05*ones(n,1),y_coor+0.05*ones(n,1),z_coor,['1';'2';'3';'4';'5';'6'],'FontSize',fs)

%% Plot trajectory
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
hold on
grid on
for ii=1:1:length(t)
jj = 6;
plot3(xx(ii,jj),yy(ii,jj),zz(ii,jj),'k.','MarkerSize',avgMs);

    if (ii==1)||(ii==200)||(ii==400)||(ii==600)||(ii==800)||(ii==length(t))
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
    hh1 = plot3(xx(1,i1),yy(1,i1),zz(1,i1),'rs','MarkerSize',1.4*avgMs);
    hh2 = plot3(xx(length(t),i1),yy(length(t),i1),zz(length(t),i1),'ro','MarkerSize',1.4*avgMs);
end
legend([hh1,hh2],'Intitial position', 'Final position')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

text(xx(length(t),:)'-0.05*ones(n,1),yy(length(t),:)'+0.05*ones(n,1),zz(length(t),:)',['1';'2';'3';'4';'5';'6'],'FontSize',fs)

%% Plot distance error %%
nn = 1:2:length(t); % Plot distance error (e) using less points

% fig one
figure
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
subplot(3,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e12(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,'Color',[0 0 1]);
plot(t(nn),e13(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,'Color',[0 0.5 0]);
plot(t(nn),e14(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,'Color',[1 0 0]);
plot(t(nn),e15(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,'Color',[0 0 0]);
plot(t(nn),e16(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,'Color',[0.75 0 0.75]);
legend('e_1_2','e_1_3','e_1_4','e_1_5','e_1_6','Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% fig two
subplot(3,1,2)

grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e23(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,'Color',[0 0 1]);
plot(t(nn),e24(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,'Color',[0 0.5 0]);
plot(t(nn),e25(nn),'-*','LineWidth',lw,'MarkerSize',1.2*avgMs,'Color',[1 0 0]);
plot(t(nn),e26(nn),'.-','LineWidth',lw,'MarkerSize',2*avgMs,'Color',[0 0 0]);
plot(t(nn),e34(nn),'-d','LineWidth',lw,'MarkerSize',0.8*avgMs,'Color',[1 0.5 0]);
legend('e_2_3','e_2_4','e_2_5','e_2_6','e_3_4','Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

% fig three
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
subplot(3,1,3)

grid on
hold on
set(gca,'Box','on','FontSize',fs)
plot(t(nn),e35(nn),'-+','LineWidth',lw,'MarkerSize',1.2*avgMs,'Color',[0 0 1]);
plot(t(nn),e36(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,'Color',[0 0.5 0]);
plot(t(nn),e45(nn),'-x','LineWidth',lw,'MarkerSize',1.4*avgMs,'Color',[0.75 0 0.75]);
plot(t(nn),e46(nn),'-s','LineWidth',lw,'MarkerSize',1.2*avgMs,'Color',[1 0 1]);
plot(t(nn),e56(nn),'-o','LineWidth',lw,'MarkerSize',0.8*avgMs,'Color',[0 0.5 0]);
legend('e_3_5','e_3_6','e_4_5','e_4_6','e_5_6','Location','EastOutside')
xlabel('Time')
ylabel('e_i_j')

%% Plot control input %%
figure
nn = 1:12:length(t);
set(gcf,'Color',[1 1 1],'Position',[232 50 560*1.8 420*1.8])
subplot(2,1,1)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),u( 1,nn),'LineWidth',1);
plot(t(nn),u( 4,nn),'LineWidth',1);
plot(t(nn),u( 7,nn),'LineWidth',1);
plot(t(nn),u(10,nn),'LineWidth',1);
plot(t(nn),u(13,nn),'LineWidth',1);
plot(t(nn),u(16,nn),'k','LineWidth',1);
xlabel('Time')
ylabel('u_i_x')

subplot(2,1,2)
grid on
hold on
set(gca,'Box','on','FontSize',0.75*fs)
plot(t(nn),u( 2,nn),'LineWidth',1);
plot(t(nn),u( 5,nn),'LineWidth',1);
plot(t(nn),u( 8,nn),'LineWidth',1);
plot(t(nn),u(11,nn),'LineWidth',1);
plot(t(nn),u(14,nn),'LineWidth',1);
plot(t(nn),u(17,nn),'k','LineWidth',1);
xlabel('Time')
ylabel('u_i_y')

% Plot errors again
figure
plot(t(nn),e12(nn),t(nn),e16(nn),t(nn),e23(nn),t(nn),e26(nn),t(nn),e34(nn),t(nn),e36(nn),t(nn),e45(nn),t(nn),e46(nn),...
     t(nn),e56(nn),t(nn),e15(nn),t(nn),e13(nn),t(nn),e14(nn),t(nn),e24(nn),t(nn),e25(nn),t(nn),e35(nn),'LineWidth',1);
xlabel('Time')
ylabel('e_i_j')
grid
xlim([0 3]);
