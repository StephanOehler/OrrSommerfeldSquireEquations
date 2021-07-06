% This file generates Figure 6 for
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

%%
% clear workspace
clear; close all
% add support files
addpath('external_files')
% options
fontsizevalue = 15;
colorlim = [-1,1] * 3.3819e-04;
% quiver options
quiv_scale = 2e3;
skipy = 10;
skipz = 5;
%%
load('LM_pregenerated/LM_data_figure6.mat')

%% create figure
% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,840, 590.625]*monitor_size(4)/1000)

% dimensioning of subplots

xl = 0.078; % left margin
xr = 0.01; % right margin
xp = (1 - xl - xr); % width of each plot
yt = 0.015; % top margin
yb = 0.13; % bottom margin
yg = [0.03,0.03,0.03,0.12]; % y gap
yp = (1 - yt - yb - sum(yg)) / (length(yg) + 1); % height of each plot
x_sub(1) = xl; % x coordinate

y_sub(5) = yb; % y coordinate
y_sub(4) = y_sub(5) + yp + yg(4);
y_sub(3) = y_sub(4) + yp + yg(3);
y_sub(2) = y_sub(3) + yp + yg(2);
y_sub(1) = y_sub(2) + yp + yg(1);

%% (a)

subplot('position', [x_sub(1) y_sub(1) xp yp]); hold on
pcolor(y,z,u_planesX.'); 
shading flat, colormap(flipud(brewermap(65,'RdBu'))) 
caxis(colorlim)

xlim([0,3*pi])
set(gca,'Box','On','linewidth',1)
set(gca, 'XTick', [0:0.5*pi:3*pi], 'XTickLabel', {'','','','','','','','','','','','','','','','',''},'Tickdir','out','TickLength',[0.0100 0.0250])
set(gca, 'YTick', [0:0.2:1], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'},'Tickdir','out','TickLength',[0.0150 0.0250])
set(gca,'FontSize',fontsizevalue)

ylabel('$z$','interpreter','latex','fontsize',fontsizevalue*sqrt(1),'rotation',0),
ylp = get(get(gca,'ylabel'), 'Position');
set(get(gca,'ylabel'), 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','Latex')
text(-pi/4,1,'(a)','interpreter','latex','fontsize',fontsizevalue)
quiver(y(1:skipy:end),z(1:skipz:end),quiv_scale * v_planesX(1:skipy:end,1:skipz:end)',quiv_scale*w_planesX(1:skipy:end,1:skipz:end)',0,'k');  hold on
ylim([0,1])
plot([0,3*pi],[0.32,0.32],'k-')

%% (b)

subplot('position', [x_sub(1) y_sub(4) xp yp]); hold on
pcolor(y,z,(uest_planesX.')); shading flat, colormap(flipud(brewermap(65,'RdBu')))
caxis(colorlim)

set(gca,'Box','On','linewidth',1)
set(gca, 'XTick', [0:0.5*pi:3*pi], 'XTickLabel', {'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$','$5\pi/2$','$3\pi$','','$4\pi$','','$5\pi$','','$6\pi$','','$7\pi$','','$8\pi$'},'Tickdir','out','TickLength',[0.0100 0.0250],'TickLabelInterpreter','Latex')
set(gca, 'YTick', [0:0.2:1], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'},'Tickdir','out','TickLength',[0.0150 0.0250])
set(gca,'FontSize',fontsizevalue)

ylabel('$z$','interpreter','latex','fontsize',fontsizevalue*sqrt(1),'rotation',0), xlabel('$y$','interpreter','latex','fontsize',fontsizevalue),
ylp = get(get(gca,'ylabel'), 'Position');
set(get(gca,'ylabel'), 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','Latex')
text(-pi/4,1,'(d)','interpreter','latex','fontsize',fontsizevalue)

quiver(y(1:skipy:end),z(1:skipz:end),quiv_scale *vest_planesX(1:skipy:end,1:skipz:end)',quiv_scale *west_planesX(1:skipy:end,1:skipz:end)',0,'k');  hold on
xlim([0,3*pi])
ylim([0,1])
plot([0,3*pi],[0.32,0.32],'k-')

%% (c)
subplot('position', [x_sub(1) y_sub(2) xp yp]); hold on
pcolor(y,z,(uest2_planesX.')); shading flat, shading flat, colormap(flipud(brewermap(65,'RdBu')))
caxis(colorlim)

set(gca,'Box','On','linewidth',1)
set(gca, 'XTick', [0:0.5*pi:3*pi], 'XTickLabel', {'','','','','','','','','','','','','','','','',''},'Tickdir','out','TickLength',[0.0100 0.0250])
set(gca, 'YTick', [0:0.2:1], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'},'Tickdir','out','TickLength',[0.0150 0.0250])
set(gca,'FontSize',fontsizevalue)

ylabel('$z$','interpreter','latex','fontsize',fontsizevalue*sqrt(1),'rotation',0),
ylp = get(get(gca,'ylabel'), 'Position');
set(get(gca,'ylabel'), 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','Latex')
text(-pi/4,1,'(b)','interpreter','latex','fontsize',fontsizevalue)

quiver(y(1:skipy:end),z(1:skipz:end),quiv_scale *vest2_planesX(1:skipy:end,1:skipz:end)',quiv_scale *west2_planesX(1:skipy:end,1:skipz:end)',0,'k');  hold on
xlim([0,3*pi])
ylim([0,1])
plot([0,3*pi],[0.32,0.32],'k-')

%% (d)
subplot('position', [x_sub(1) y_sub(3) xp yp]); hold on
pcolor(y,z,(uest3_planesX.')); shading flat, shading flat, colormap(flipud(brewermap(65,'RdBu'))) % plot intial frame
caxis(colorlim)

set(gca,'Box','On','linewidth',1)
set(gca, 'XTick', [0:0.5*pi:3*pi], 'XTickLabel', {'','','','','','','','','','','','','','','','',''},'Tickdir','out','TickLength',[0.0100 0.0250])
set(gca, 'YTick', [0:0.2:1], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'},'Tickdir','out','TickLength',[0.0150 0.0250])
set(gca,'FontSize',fontsizevalue)
ylb = ylabel('$z$','interpreter','latex','fontsize',fontsizevalue*sqrt(1)); 
set(ylb, 'Rotation',0,'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','Latex')

text(-pi/4,1,'(c)','interpreter','latex','fontsize',fontsizevalue)

quiver(y(1:skipy:end),z(1:skipz:end),quiv_scale *vest3_planesX(1:skipy:end,1:skipz:end)',quiv_scale *west3_planesX(1:skipy:end,1:skipz:end)',0,'k');  hold on
ylim([0,1])
xlim([0,3*pi])
plot([0,3*pi],[0.32,0.32],'k-')

%% save figure (non-vectorized)
set(gcf,'renderer','opengl')

