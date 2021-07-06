% This file generates Figure 5 for 
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735


% clear workspace
clear; close all
% add support files
addpath('external_files')

%% setup figure

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,2*sqrt(2),1]*monitor_size(4)/2)

% scale fonts and lines
textfontsize = 14;
fontsizevalue = 16;

% define axis limit
caxislimits = [0.2,0.9];


fdim13.xl = 0.06; % left margin
fdim13.xr = 0.005; % right margin
fdim13.xg = [0.012,0.012];
fdim13.xp = (1 - fdim13.xl - fdim13.xr - sum(fdim13.xg))/(length(fdim13.xg) + 1); % width of each plot
fdim13.yt = 0.09; % top margin
fdim13.yb = 0.15; % bottom margin
fdim13.yp = 1 - fdim13.yb - fdim13.yt; % width of each plot
fdim13.x(1) = fdim13.xl; % x coordinate of (a)
fdim13.x(2) = fdim13.x(1) + fdim13.xg(1) + fdim13.xp; % x coordinate of (b)
fdim13.x(3) = fdim13.x(2) + fdim13.xg(2) + fdim13.xp; % x coordinate of (c)
fdim13.y(1) = fdim13.yb; % y coordinate



load('fig5/fig5.mat')

% plot figure (c)
subplot('position', [fdim13.x(3) fdim13.y(1) fdim13.xp fdim13.yp]), hold on
contourf(KY,KX,(real(IOC_gam_z./P_gam_mat))); view(90,-90)
axis tight
set(gca,'TickLabelInterpreter','Latex')
set(gca,'FontSize',textfontsize)
ylabel('$k_x$','interpreter','latex')
set(gca,'YTick',[0.0:1:8])
set(gca,'XTick',[0.0:2:21],'XTickLabel',{})
set(gca,'TickDir','Out','TickLength',[0.0200 0.0250],'box','on')
text(22.2,0.1,'(c)','interpreter','latex','interpreter','latex','fontsize',fontsizevalue)
colormap(brewermap(75,'YlGnBu'))
caxis([0.24,0.98])
plot([0,8],[0,8],'w--','linewidth',0.5) % plot line of symmetry
shading interp
set(get(gca,'Title'),'Position',[22,4,0])

% plot figure (a)
subplot('position', [fdim13.x(1) fdim13.y(1) fdim13.xp fdim13.yp]), hold on
contourf(KY,KX,(real(OE_gam_u./P_gam_mat))); view(90,-90)
axis tight
set(gca,'TickLabelInterpreter','Latex')
set(gca,'FontSize',textfontsize)
ylabel('$k_x$','interpreter','latex')
xlabel('$k_y$','interpreter','latex')
set(gca,'YTick',[0.0:1:8])
set(gca,'XTick',[0.0:2:21])
set(gca,'TickDir','Out','TickLength',[0.0200 0.0250],'box','on')
text(22.2,0.1,'(a)','interpreter','latex','interpreter','latex','fontsize',fontsizevalue)
set(get(gca,'Title'),'Position',[22,4,0])

% plot line of symmetry
colormap(brewermap(75,'YlGnBu'))
caxis([0.24,0.98])
plot([0,8],[0,8],'w--','linewidth',0.5)
shading interp


% plot figure (b)
subplot('position', [fdim13.x(2) fdim13.y(1) fdim13.xp fdim13.yp]), hold on
contourf(KY,KX,(real(FIC_gam_z./P_gam_mat))); view(90,-90)
axis tight
set(gca,'TickLabelInterpreter','Latex')
set(gca,'FontSize',textfontsize)
ylabel('$k_x$','interpreter','latex')
set(gca,'YTick',[0.0:1:8])
set(gca,'XTick',[0.0:2:21],'XTickLabel',{})
set(gca,'TickDir','Out','TickLength',[0.0200 0.0250],'box','on')
text(22.2,0.1,'(b)','interpreter','latex','interpreter','latex','fontsize',fontsizevalue)
set(get(gca,'Title'),'Position',[22,4,0])

% plot line of symmetry
colormap(brewermap(75,'YlGnBu'))
caxis([0.24,0.98])
plot([0,8],[0,8],'w--','linewidth',0.5)
shading interp

% plot colorbar
c = colorbar;
c.Location = 'North';
c.Ticks = [0.2:0.1:1];
c.Color = [1,1,1];

% switch to non-vector renderer
set(gcf,'renderer','opengl')

% save figure
print('fig5','-depsc','-r200')

