% This file generates Figure 2 for 
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735


% clear workspace
clear; close all
% add support files
addpath('external_files')

%% setup figure

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,2,1.1*1/sqrt(2)]*monitor_size(4)/1.2)

% scale fonts and lines
textfontsize = 14;
linefont = 1.2;
lightgreen = [0.2745,0.7373,1];%[0 0.4470 0.7410];

% define axis limit
caxislimits = [0.2,0.9];

Nx = 200; % select number of chebyshev grid points to be selected
Nboth = 200;
Re = 2000;
%%

%flowDNS = chanDNS.Madrid2000(1,1); % form linear model
flow = channelOSS.StateSpace(1,1,Re,Nx,'turbulent','top',Nboth);

iz_uvw = [1:Nboth+1;Nboth+2:2*Nboth+2;2*Nboth+3:3*Nboth+3];
iz_u = 1:Nboth+1;
iz_v = Nboth+2:2*Nboth+2;
iz_w = 2*Nboth+3:3*Nboth+3;



%%
fdim12v2.xl = 0.06; % left margin
fdim12v2.xr = 0.01; % right margin
fdim12v2.xg = [0.12]; % x-gap
fdim12v2.xp = (1 - fdim12v2.xl - fdim12v2.xr - sum(fdim12v2.xg))/(length(fdim12v2.xg) + 1); % width of each plot
fdim12v2.yt = 0.14; % top margin
fdim12v2.yb = 0.135; % bottom margin
fdim12v2.yp = 1 - fdim12v2.yb - fdim12v2.yt; % width of each plot
% x and y dimension of plot
fdim12v2.x(1) = fdim12v2.xl;
fdim12v2.x(2) = fdim12v2.x(1) + fdim12v2.xg(1) + fdim12v2.xp;
fdim12v2.y(1) = fdim12v2.yb;

%%


%% load in DNS data (uncontrolled)
load('DNS_pregenerated/DNS_norms')

P2_DNS_rms = (flow.w_out .\ (sum(P2_DNSuvw(iz_uvw),1).'));
P2_DNS_rms_u = ((flow.w_out .\ P2_DNSuvw(iz_u)));
P2_DNS_rms_v = ((flow.w_out .\ P2_DNSuvw(iz_v)));
P2_DNS_rms_w = ((flow.w_out .\ P2_DNSuvw(iz_w)));

P2_DNS_max = max(P2_DNS_rms);

%% load in LM data (uncontrolled)
load('fig8/fig8.mat','P_rms_output_mat')

P_rms = squeeze(sum(P_rms_output_mat,[1 2]));

P2_sim_rms = (flow.w_out .\ (sum(P_rms(iz_uvw),1).'));
P2_sim_rms_u = ((flow.w_out .\ P_rms(iz_u)));
P2_sim_rms_v = ((flow.w_out .\ P_rms(iz_v)));
P2_sim_rms_w = ((flow.w_out .\ P_rms(iz_w)));
P2_max = max(P2_sim_rms);

%% plot lines
subplot('position', [fdim12v2.x(2) fdim12v2.y(1) fdim12v2.xp fdim12v2.yp]), hold on
hold on

% DNS
plot(1-flow.z_out,P2_DNS_rms/P2_DNS_max,'-','color',lightgreen,'linewidth',linefont)
plot(1-flow.z_out,P2_DNS_rms_u/P2_DNS_max,'--','color',lightgreen,'linewidth',linefont)
plot(1-flow.z_out,P2_DNS_rms_v/P2_DNS_max,'-.','color',lightgreen,'linewidth',linefont)
plot(1-flow.z_out,P2_DNS_rms_w/P2_DNS_max,':','color',lightgreen,'linewidth',linefont)

% LM
plot(1-flow.z_out,P2_sim_rms/P2_max,'k-','linewidth',linefont)
plot(1-flow.z_out,P2_sim_rms_u/P2_max,'k--','linewidth',linefont)
plot(1-flow.z_out,P2_sim_rms_v/P2_max,'k-.','linewidth',linefont)
plot(1-flow.z_out,P2_sim_rms_w/P2_max,'k:','linewidth',linefont)



% annotate:
xlabel('$z$','interpreter','latex','FontSize',textfontsize), 
set(gca,'FontSize',textfontsize)
set(gca,'YTick',[0:0.2:1])
set(gca,'XTick',[0:0.2:1])
axis tight 
ylim([0,1])
set(gca,'TickLabelInterpreter','Latex')
box on
set(gca,'TickLabelInterpreter','latex')

% apply custom ticks
TicksforLabel = {'$0$','','','','','$\frac{1}{4}\max\|\hat{\textit{\textbf u}}\|_2^2$','','','','','$\frac{1}{2}\max\|\hat{\textit{\textbf u}}\|_2^2$','','','','','$\frac{3}{4}\max\|\hat{\textit{\textbf u}}\|_2^2$','','','','','$\max\|\hat{\textit{\textbf u}}\|_2^2$'};
set(gca,'YTick',[0.0:0.05:1],'YTickLabel',TicksforLabel,'Tickdir','out','TickLength',[0.0100 0.0250])
TicksforLabel = {'$0$','','$0.1$','','$0.2$','','$0.3$','','$0.4$','','$0.5$','','$0.6$','','$0.7$','','$0.8$','','$0.9$','','$1$'};
set(gca,'XTick',[0.0:0.05:1],'XTickLabel',TicksforLabel,'Tickdir','out','TickLength',[0.0100 0.0250])

% save as vector graphics
set(gcf,'renderer','painters')
print('-depsc','fig2','-r500')

