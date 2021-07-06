% This file generates Figure 8 for
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735


% clear workspace
clear; close all
% add support files
addpath('external_files')

%% setup figure

% Scale figure
monitor_size = get(0,'screensize');
figure('color','white','position',[0,0,2,1]*monitor_size(4)/2)


% scale fonts and lines
textfontsize = 16;
linefont = 1.2;

% define axis limit
caxislimits = [0.2,0.9];

% Change tick labelsS
TicksforLabelx = {'','0.1','','0.3','','0.5','','0.7','','0.9',''};
TicksforLabely = {'','90\%','','70\%','','50\%','','30\%','','10\%',''};
emptyTicksforLabelx = {'','','','','','','','','','',''};
emptyTicksforLabely = {'','','','','','','','','','',''};

% place figures
fdim22.xl = 0.07; % left margin
fdim22.xr = 0.09; % right margin
fdim22.xg = [0.005]; % x-gap
fdim22.xp = (1 - fdim22.xl - fdim22.xr - sum(fdim22.xg))/(length(fdim22.xg) + 1); % width of each plot
fdim22.yg = [0.01]; % y-gap
fdim22.yt = 0.06; % top margin
fdim22.yb = 0.11; % bottom margin
%yp = (1 - yt - yb - sum(yg))/(length(yg) + 1); % height of each plot
fdim22.yp = (1 - fdim22.yt - fdim22.yb - sum(fdim22.yg)) / ( length(fdim22.yg) + 1 ); % height of each plot

% x-coordinate
fdim22.x(1) = fdim22.xl;
fdim22.x(2) = fdim22.x(1) + fdim22.xg(1) + fdim22.xp;
fdim22.x(3) = fdim22.xl;
fdim22.x(4) = fdim22.x(1) + fdim22.xg(1) + fdim22.xp;

% y-coordinate
fdim22.y(3) = fdim22.yb;
fdim22.y(1) = fdim22.y(3) + fdim22.yp + fdim22.yg(1);
fdim22.y(4) = fdim22.yb;
fdim22.y(2) = fdim22.y(3) + fdim22.yp + fdim22.yg(1);



% load data
load('fig8/fig8.mat')

flow = channelOSS.StateSpace(1,1,Re,Nx,'turbulent','top',Nboth);
zgridhalf = flow.z_out;

% select displayed data: either all measurements, u, v, or w
var = {1:3,1,2,3}; %{'uvw','u','v','w'};
% indicies corresponding to u v and w (1st,2nd and 3rd row)
idx_iii =  [1:Nboth+1;Nboth+2:2*Nboth+2;2*Nboth+3:3*Nboth+3];


% loop over the four cases: {'uvw','u','v','w'};
for ijk = 1:4
    
    if ijk == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp = squeeze(sum(sum(P_rms_output_mat([1:end,2:end],[1:end,2:end],:),2),1));
        P_rms = sum(temp(idx_iii(var{ijk},:)),1);
        P_norm = sum(temp(idx_iii(var{ijk},:)) ,1) ./ flow.w_out.';
        temp = squeeze(sum(sum(OE_rms_u([1:end,2:end],[1:end,2:end],:),2),1)) ;
        rms_u = sum(temp(idx_iii(var{ijk},:)),1) ./ P_rms;
        temp = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],:),2),1)) ;
        rms_z = sum(temp(idx_iii(var{ijk},:)),1) ./ P_rms;
        temp = squeeze(sum(sum(IOC_rms_uz([1:end,2:end],[1:end,2:end],:),2),1)) ;
        rms_uz = sum(temp(idx_iii(var{ijk},:)),1) ./ P_rms;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp = squeeze(sum(sum(P_rms_output_mat([1:end,2:end],[1:end,2:end],:),2),1));
        P_rms = (temp(idx_iii(var{ijk},:)));
        temp = squeeze(sum(sum(OE_rms_u([1:end,2:end],[1:end,2:end],:),2),1)) ;
        rms_u = (temp(idx_iii(var{ijk},:))) ./ P_rms;
        temp = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],:),2),1)) ;
        rms_z = (temp(idx_iii(var{ijk},:))) ./ P_rms;
        temp = squeeze(sum(sum(IOC_rms_uz([1:end,2:end],[1:end,2:end],:),2),1)) ;
        rms_uz = (temp(idx_iii(var{ijk},:))) ./ P_rms;
        P_norm = flow.w_out .\ P_rms;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % create subplot
    subplot('position', [fdim22.x(ijk) fdim22.y(ijk) fdim22.xp fdim22.yp]), hold on
    
    % create blue plot with y-axis on the right
    yyaxis right
    % plot blue data
    plot(1-zgridhalf,(P_norm)/max(P_norm),'-','color',[0.2745,0.7373,1.0000],'linewidth',linefont),
    set(gca,'YColor',[0.2745,0.7373,1.0000])
    
    % change ticks
    set(gca,'Tickdir','out','TickLength',[0.0100 0.0250])
   
    % create ylabel
    if any(ijk == [2,4])
        ylabel('$\mathbf{E}_z$','interpreter','latex','FontSize',textfontsize)
        set(gca,'YTick',[0.0:0.1:1],'YTickLabel',TicksforLabelx)
    else
        set(gca,'YTick',[0.0:0.1:1],'YTickLabel',emptyTicksforLabelx)
    end
    
    % give each subplot a leter
    if ijk == 1
        text(0.93, 0.93,{'(a)'},'interpreter','latex','FontSize',textfontsize)
    elseif ijk == 2
        text(0.93, 0.93,{'(b)'},'interpreter','latex','FontSize',textfontsize)
    elseif ijk == 3
        text(0.93, 0.07,{'(c)'},'interpreter','latex','FontSize',textfontsize)
    elseif ijk == 4
        text(0.93, 0.07,{'(d)'},'interpreter','latex','FontSize',textfontsize)
    end
    
    % plot main data (black)
    yyaxis left
    hold on, plot(1-zgridhalf,real(rms_u),'--','color','k','linewidth',linefont),
    hold on, plot(1-zgridhalf,real(rms_uz),'-','color','k','linewidth',linefont),
    hold on, plot(1-zgridhalf,real(rms_z),'-.','color','k','linewidth',linefont),
    
    
    % style figure
    set(gca,'FontSize',textfontsize,'Ycolor',[0,0,0])
    ylim([0.0,1.0])
    set(gca,'TickLabelInterpreter','Latex')
    box on
    
    % change Xticks
    if any(ijk == [3,4])
        xlabel('$z$','interpreter','latex','FontSize',textfontsize)
        set(gca,'XTick',[0.0:0.1:1],'XTickLabel',TicksforLabelx)
    else
        set(gca,'XTick',[0.0:0.1:1],'XTickLabel',emptyTicksforLabelx)
    end
    
    % change Yticks
    if any(ijk == [1,3])
        ylabel('$\epsilon$','interpreter','latex','FontSize',textfontsize) %: reduction of $\mathbf{E}_K$
        set(gca,'YTick',[0.0:0.1:1],'YTickLabel',TicksforLabely)
    else
        set(gca,'YTick',[0.0:0.1:1],'YTickLabel',emptyTicksforLabely)
    end
      
    
end

% save figure
print('-depsc2','fig8')

