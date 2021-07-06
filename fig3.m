% This file generates Figure 3 for
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

%%
% clear workspace
clear; close all
% add support files
addpath('external_files')
% options
fontsizevalue = 14;

%% create figure
% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,425,446.25]*monitor_size(4)/800)


% dimensioning of subplots

xl = 0.118; % left margin
xr = 0.04; % right margin
xp = (1 - xl - xr); % width of each plot
yt = 0.015; % top margin
yb = 0.11; % bottom margin
yg = [0.03,0.03,0.12];
yp = (1 - yt - yb - sum(yg)) / (length(yg) + 1); % height of each plot
x_sub(1) = xl;

y_sub(4) = yb;
y_sub(3) = y_sub(4) + yp + yg(3);
y_sub(2) = y_sub(3) + yp + yg(2);
y_sub(1) = y_sub(2) + yp + yg(1);

%% plot DNS
load('DNS_pregenerated/DNS_data_figure3.mat')

subplot('position', [x_sub(1) y_sub(2) xp yp]); hold on
surf(y(1:end/2),z,-u_planesX(1:end/2,:)'); shading flat, colormap(brewermap(128,'RdBu'))
colorlim = [-2.5,2.5];
caxis(colorlim)
xlim([0,1.5*pi])

set(gca,'Box','On','linewidth',1)
set(gca, 'YTick', [0:0.2:1], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'},'Tickdir','out','TickLength',[0.0150 0.0250])
set(gca,'FontSize',fontsizevalue)

set(gca, 'XTick', [0:0.5*pi:3*pi], 'XTickLabel', {'0','0.5$\pi$','$\pi$','1.5$\pi$','2$\pi$','2.5$\pi$','3$\pi$','','4\pi','','5\pi','','6\pi','','7\pi','','8\pi'},'Tickdir','out','TickLength',[0.0100 0.0250])
xlabel('$y$','interpreter','latex','fontsize',fontsizevalue),

ylabel('$z$','interpreter','latex','fontsize',fontsizevalue,'rotation',0),
ylp = get(get(gca,'ylabel'), 'Position');
set(get(gca,'ylabel'), 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','Latex')
text(-.65,1,'(b)','interpreter','latex')
%% LM
load('LM_pregenerated/LM_data_figure3')

subplot('position', [x_sub(1) y_sub(1) xp yp]); hold on
surf(y(1:end/2),z,-u_planesX(1:end/2,:)'); shading flat, colormap(brewermap(128,'RdBu'))
colorlim = [-2.5,2.5]*1e-4;
caxis(colorlim)

xlim([0,1.5*pi])
set(gca,'Box','On','linewidth',1)
set(gca, 'XTick', [0:0.5*pi:3*pi], 'XTickLabel', {'','','','','','','','','','','','','','','','',''},'Tickdir','out','TickLength',[0.0100 0.0250])
set(gca, 'YTick', [0:0.2:1], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'},'Tickdir','out','TickLength',[0.0150 0.0250])
set(gca,'FontSize',fontsizevalue)

ylabel('$z$','interpreter','latex','fontsize',fontsizevalue,'rotation',0),
ylp = get(get(gca,'ylabel'), 'Position');
set(get(gca,'ylabel'), 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','Latex')
text(-.65,1,'(a)','interpreter','latex')

%% save figure (non-vectorized)
set(gcf,'renderer','opengl')
print('-depsc','fig3','-r1200')





