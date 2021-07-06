% This file generates Figure 1 for 
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


% LM data
load('fig5/fig5.mat')

%%

fdim13.xl = 0.06; % left margin
fdim13.xr = 0.005; % right margin
fdim13.xg = [0.012,0.012];
fdim13.xp = (1 - fdim13.xl - fdim13.xr - sum(fdim13.xg))/(length(fdim13.xg) + 1); % width of each plot
fdim13.yt = 0.09; % top margin
fdim13.yb = 0.15; % bottom margin
fdim13.yp = 1 - fdim13.yb - fdim13.yt; % width of each plot
fdim13.x(1) = fdim13.xl; % x coordinate of (a)
fdim13.x(2) = fdim13.x(1) + fdim13.xg(1) + fdim13.xp; % x coordinate of (b)
fdim13.x(3) = fdim13.x(2) + fdim13.xg(2) + fdim13.xp; % x coordinate of the colourbar
fdim13.y(1) = fdim13.yb; % y coordinate

% contour range
range = -2.4:0.2:0;
caxisvalues = [-2.4,0];

%% plot LM data
subplot('position', [fdim13.x(1) fdim13.y(1) fdim13.xp fdim13.yp]), hold on
contourf(KY,KX,log10(real(P_gam_mat)/max(max(real(P_gam_mat)))),range); view(90,-90)

axis tight
set(gca,'TickLabelInterpreter','Latex')
set(gca,'FontSize',textfontsize)
ylabel('$k_x$','interpreter','latex')
xlabel('$k_y$','interpreter','latex')
set(gca,'YTick',[0.0:1:8])
set(gca,'XTick',[0.0:2:21])
set(gca,'TickDir','Out','TickLength',[0.0150 0.0250])

caxis(caxisvalues)
colormap((brewermap(59,'YlGnBu')))
plot([0,8],[0,8],'k--','linewidth',0.5)
plot([0,6+2/3/5,6+2/3/5,0,0],[0,0,0.5+0.25/5,0.5+0.25/5,0],'w','linewidth',0.5);
shading interp
text(22.2,3.6,'$\textbf{LM}$','interpreter','latex','interpreter','latex','fontsize',fontsizevalue)


%% load in DNS results and then take the average average

load DNS_pregenerated\P2kxky_DNS_large.mat
P2average = P2kxky_DNSuvw(1:33,1:33);
P2average(2:33,1:33) = (P2kxky_DNSuvw(2:33,1:33) + P2kxky_DNSuvw(65:-1:34,1:33))/2;
P2average(1,1) = NaN;
%% Plot DNS results

subplot('position', [fdim13.x(2) fdim13.y(1) fdim13.xp fdim13.yp]), hold on
contourf(KY,KX,log10(real(P2average)/max(max(real(P2average)))),range); view(90,-90)

axis tight
set(gca,'TickLabelInterpreter','Latex')
set(gca,'FontSize',textfontsize)
ylabel('$k_x$','interpreter','latex')
set(gca,'YTick',[0.0:1:8])
set(gca,'XTick',[0.0:2:21])
set(gca,'XTickLabel',{'','','','','','','','','','',''})
set(gca,'TickDir','Out','TickLength',[0.0150 0.0250])

caxis(caxisvalues)
colormap((brewermap(59,'YlGnBu')))
plot([0,8],[0,8],'k--','linewidth',0.5)
plot([0,6+2/3/5,6+2/3/5,0,0],[0,0,0.5+0.25/5,0.5+0.25/5,0],'w','linewidth',0.5);
shading interp
text(22.2,3.3,'$\textbf{DNS}$','interpreter','latex','interpreter','latex','fontsize',fontsizevalue)


%% plot colourbar
c2 = colorbar;
c2.Position =  [fdim13.x(3) fdim13.y(1) 0.02 fdim13.yp];
c2.Ticks = range;
c2.TickLabelInterpreter = 'latex';
c2.TickLabels = {'','$10^{-2.2}$','','$10^{-1.8}$','','$10^{-1.4}$','','$10^{-1}$','','$10^{-0.6}$','','$10^{-0.2}$',''};


%% save figure

set(gcf,'renderer','opengl')
print('fig1','-depsc','-r1200')


