% This file generates the support files for Figure 5 of
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

clear

addpath('..')
%% User options


% select number of chebyshev grid points to be selected
Nx = 200; % Number of grid points of the state space model
Nboth = 200; % Number of grid points of the output

% select part of channel that is to be estimated
KX = 0:0.25:8; % range of kx
KY = 0:2/3:64/3; % range of ky

% Reynolds number of the DNS data
Re = 2000; % Re_\tau

% generates a State space model of the turbulent channel flow
flow = channelOSS.StateSpace(1,1,Re,Nx,'turbulent','top',Nboth);

%% Setup model and sensors / estimation locations

% control cost
Rp_vel = (eye(3) * 1e-4)^2;
Rp_wall = (eye(2) * 1e-4)^2;

% sensor cost
Vp_vel = (eye(3) * 1e-4)^2;
Vwall = (eye(2) * 1e-4)^2;

%% simulation results at each wavenumber pair

% energy of system
P_gam_mat = zeros(length(KX),length(KY));

% Full information control (point actuator / blowing and suction)
FIC_gam_z = zeros(length(KX),length(KY));
FIC_gam_BBs = zeros(length(KX),length(KY));

% Optimal Estimation / Actuation everywhere (point sensor / shear stress at wall)
OE_gam_u = zeros(length(KX),length(KY));
OE_gam_taux = zeros(length(KX),length(KY));

% Input output control (in flow / at wall)
IOC_gam_z = zeros(length(KX),length(KY));
IOC_gam_bbs = zeros(length(KX),length(KY));

%% OE result for sensor placement for optimal blow suction IOC with in-flow sensor

% sensor and actuator gaussian width
var_a = 0.02;
var_s = 0.02;

meas_loc_s = 0.68;% 1 - 0.32
meas_loc_a = 0.68;% 1 - 0.32

% single wall height actuator
Bint_dz = exp( -((flow.geom.x - meas_loc_a) / var_a ) .^2 );
Bint_dz = [Bint_dz, zeros(size(Bint_dz)), zeros(size(Bint_dz)); ...
     zeros(size(Bint_dz)), Bint_dz, zeros(size(Bint_dz)); ...
     zeros(size(Bint_dz)), zeros(size(Bint_dz)), Bint_dz];
 
% single wall height sensor 
Cint_u = ( diag(flow.w_out) * exp( -((flow.z_out - meas_loc_s) / var_s ) .^2 ) ).';
Cint_u = [Cint_u, zeros(size(Cint_u)), zeros(size(Cint_u)); ...
     zeros(size(Cint_u)), Cint_u, zeros(size(Cint_u)); ...
     zeros(size(Cint_u)), zeros(size(Cint_u)), Cint_u];

%% loop over kx and ky
for i = 1:length(KX)
    tic
    kx = KX(i); % loop over ky
    parfor j = 1:length(KY)
        tic
        % linear model at this wavenumber pair
        flow = channelOSS.StateSpace(kx,KY(j),Re,Nx,'turbulent','top',Nboth);
        C1 = flow.Cw;
        %%
        Z1 = lyap(flow.A', C1' * C1);
        P_gam_mat(i,j) = real(trace(flow.Bw' * Z1 * flow.Bw));
        %%
        C2 = Cint_u * flow.C;
        [Y1,~,~] = care(flow.A',C2',flow.Bw*flow.Bw',Vp_vel);
        OE_gam_u(i,j) = trace(C1 * Y1 * C1');
        %%
        C2_wall = flow.Ctaup(1:2,:);
        [Y2,~,~] = care(flow.A',C2_wall',flow.Bw*flow.Bw',Vwall);
        OE_gam_taux(i,j) = trace(C1 * Y2 * C1');
        %%
        B2_1 = flow.B * Bint_dz;      
        [X1,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel); 
        FIC_gam_z(i,j) = real(trace(flow.Bw' * X1 * flow.Bw));
        %%
        B2_2 = flow.Bbs(:,[1,7]);
        [X2,~,~] = care(flow.A,B2_2,C1'*C1,Rp_wall);
        FIC_gam_BBs(i,j) = real(trace(flow.Bw' * X2 * flow.Bw));
        %%      
        IOC_gam_z(i,j) = trace(flow.Bw' * X1 * flow.Bw) + trace( inv(Rp_vel) * B2_1' * X1*Y1*X1 * B2_1);
        IOC_gam_bbs(i,j) = trace(flow.Bw' * X2 * flow.Bw) + trace( inv(Rp_wall) * B2_2' * X2*Y2*X2 * B2_2);

        disp(['i = ',num2str(i),', j = ',num2str(j),', toc = ',num2str(toc)])
    end
    toc
    save('fig5')
end
%% overall energgy of the flow
P_gam = sum(sum(P_gam_mat([1:end,2:end],[1:end,2:end]),2),1);

%% overal energy norms

gam_u = sum(sum(OE_gam_u([1:end,2:end],[1:end,2:end]),2),1) / P_gam;

gam_taux = sum(sum(OE_gam_taux([1:end,2:end],[1:end,2:end]),2),1) / P_gam;

gam_z = sum(sum(FIC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;

gam_bbs = sum(sum(FIC_gam_BBs([1:end,2:end],[1:end,2:end]),2),1) / P_gam;

gam_uz = sum(sum(IOC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;

gam_ubbs = sum(sum(IOC_gam_bbs([1:end,2:end],[1:end,2:end]),2),1) / P_gam;

%% save data
save('fig5','gam_u','gam_taux','gam_z','gam_bbs','gam_uz','gam_ubbs', ...
    'P_gam_mat','FIC_gam_z','FIC_gam_BBs','OE_gam_u','OE_gam_taux', ...
    'IOC_gam_z','IOC_gam_bbs','Re','Nx','Nboth','KX','KY',...
    'var_a','var_s','meas_loc_s','meas_loc_a')

