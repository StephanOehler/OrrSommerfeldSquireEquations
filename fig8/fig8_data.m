% This file generates the support files for Figure 8 of
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

clear

addpath('..')
%% User options


% select number of chebyshev grid points to be selected
Nx = 200; % Number of grid points of the state space model
Nboth = 200; % Number of grid points of the output

% select part of channel that is to be estimated
KX = 0:0.25:0.5; % range of kx
KY = 0:2/3:6; % range of ky


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

%% allocate memory

% z--> in-flow actuator
% u--> in-flow sensor
% bbs--> blowing and suction actuator
% tau--> shear stress sensor (streamwise shear)

% Allocate memory for uncontrolled flow
P_rms_output_mat = zeros(length(KX),length(KY),  3*(flow.N_out + 1));
P_gam_mat = zeros(length(KX),length(KY));

% RMS for FIC, single plane actuator and blowing and suction
FIC_rms_z = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + 3);
FIC_rms_bbs = zeros(length(KX),length(KY), 3*(flow.N_out + 1) + 2);

% H2 norm for FIC, single plane actuator and blowing and suction
FIC_gam_z = zeros(length(KX),length(KY));
FIC_gam_bbs = zeros(length(KX),length(KY));

% Repeat for IOC
IOC_gam_uz = zeros(length(KX),length(KY));
IOC_rms_uz = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + 3);
IOC_gam_taubbs = zeros(length(KX),length(KY));
IOC_rms_taubbs = zeros(length(KX),length(KY), 3*(flow.N_out + 1) + 2);

% Repeat for Optimal Estimation / Actuation Everywhere
OE_gam_u = zeros(length(KX),length(KY));
OE_rms_u = zeros(length(KX),length(KY),  3*(flow.N_out + 1));
AE_rms_forcing = zeros(length(KX),length(KY),  3*(flow.N_out + 1));



OE_gam_tau = zeros(length(KX),length(KY));
OE_rms_tau = zeros(length(KX),length(KY),  3*(flow.N_out + 1));

% IOC result for optimal in-flow sensor placement for blowing suction IOC
IOC_gam_ubbs = zeros(length(KX),length(KY));
IOC_rms_ubbs = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + 3);

%%
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
    kx = KX(i);
    parfor j = 1:length(KY)
        tic
        % linear model at this wavenumber pair
        flow = channelBC.StateSpace(kx,KY(j),Re,Nx,'turbulent','top',Nboth);
        C1 = flow.Cw;
        
        %% uncontrolled
        Z1 = lyap(flow.A', C1' * C1);
        P_gam_mat(i,j) = real(trace(flow.Bw' * Z1 * flow.Bw));
        Z2 = lyap(flow.A, flow.Bw * flow.Bw');
        P_rms_output_mat(i,j,:) = real(diag(C1 * Z2 * C1'));
        %% OE in-flow
        C2 = Cint_u * flow.C;
        [Y,~,~] = care(flow.A',C2',flow.Bw*flow.Bw',Vp_vel);
        OE_gam_u(i,j) = trace(C1 * Y * C1');
        OE_rms_u(i,j,:) = diag(C1 * Y * C1'); 
            
        % forcing energy for ideal AE case
        C3 = C1 * (Y * C2' / Vp_vel) * C2; % Cw * L * Cy
        AE_rms_forcing(i,j,:) = diag(C3 * Y * C3');        
        %% OE at-wall
        C2 = flow.Ctaup(1:2,:);
        [Y,~,~] = care(flow.A',C2',flow.Bw*flow.Bw',Vwall);
        OE_gam_tau(i,j) = trace(C1 * Y * C1');
        OE_rms_tau(i,j,:) = diag(C1 * Y * C1');
        
        %% FIC in-flow
        B2_1 = flow.B * Bint_dz;      
        [X1,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel);
        FIC_gam_z(i,j) = real(trace(flow.Bw' * X1 * flow.Bw));
        
        F = Rp_vel \ B2_1' * X1;
        Z_a = flow.A - B2_1 * F;
        Z_b = flow.Bw;
        Z_c = [C1;-sqrt(Rp_vel) * F];
        X2 = lyap(Z_a, Z_b * Z_b');
        FIC_rms_z(i,j,:) = real(diag(Z_c * X2 * Z_c'));
        
        %% FIC at-wall
        B2_2 = flow.Bbs(:,[1,7]);
        [X3,~,~] = care(flow.A,B2_2,C1'*C1,Rp_wall);
        FIC_gam_bbs(i,j) = real(trace(flow.Bw' * X3 * flow.Bw));
        
        F = Rp_wall \ B2_2' * X3;
        Z_a = flow.A - B2_2 * F;
        Z_b = flow.Bw;
        Z_c = [C1;-sqrt(Rp_wall) * F];
        X4 = lyap(Z_a, Z_b * Z_b');
        FIC_rms_bbs(i,j,:) = real(diag(Z_c * X4 * Z_c')); 
        
        %% IOC in-flow
        B2_1 = flow.B * Bint_dz;
        [X,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel);
        F = Rp_vel \ B2_1' * X;
        C2 = Cint_u * flow.C;
        [Y,~,~] = care(flow.A',C2',flow.Bw*flow.Bw',Vp_vel);
        L = Y * C2' / Vp_vel;
        IOC_gam_uz(i,j) = trace(flow.Bw' * X * flow.Bw) + trace( inv(Rp_vel) * B2_1' * X*Y*X * B2_1);
               
        Z_a = [flow.A, -B2_1*F; L * C2, flow.A - B2_1 * F - L * C2];
        Z_b = [flow.Bw,zeros(size(flow.Bw,1),size(L*Vp_vel,2)) ; zeros(size(L*Vp_vel,1),size(flow.Bw,2)), L * sqrt(Vp_vel)]; % , zeros(size(flow.Bw,1),size(L*Vp_vel,2))];%; zeros(size(L*Vp05,1),size(flow.Bw,2)), L * Vp05];
        Z_c = [C1, zeros(size(C1,1),size(Rp_vel*F,2)); zeros(size(Rp_vel*F,1),size(C1,2)), -sqrt(Rp_vel) * F];
        Z = lyap(Z_a, Z_b * Z_b');
        IOC_rms_uz(i,j,:) = real(diag(Z_c * Z * Z_c'));
        
        %% IOC mixed
        B2_1 = flow.B * Bint_dz;
        [X,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel);
        F = Rp_vel \ B2_1' * X;
        C2_2 = flow.Ctaup(1:2,:);
        [Y,~,~] = care(flow.A',C2_2',flow.Bw*flow.Bw',Vwall);
        L = Y * C2_2' / Vwall;
        IOC_gam_ubbs(i,j) = trace(flow.Bw' * X * flow.Bw) + trace( inv(Rp_vel) * B2_1' * X*Y*X * B2_1);
               
        Z_a = [flow.A, -B2_1*F; L * C2_2, flow.A - B2_1 * F - L * C2_2];
        Z_b = [flow.Bw,zeros(size(flow.Bw,1),size(L*Vwall,2)) ; zeros(size(L*Vwall,1),size(flow.Bw,2)), L * sqrt(Vwall)]; % , zeros(size(flow.Bw,1),size(L*Vp_vel,2))];%; zeros(size(L*Vp05,1),size(flow.Bw,2)), L * Vp05];
        Z_c = [C1, zeros(size(C1,1),size(Rp_vel*F,2)); zeros(size(Rp_vel*F,1),size(C1,2)), -sqrt(Rp_vel) * F];
        Z = lyap(Z_a, Z_b * Z_b');
        IOC_rms_ubbs(i,j,:) = real(diag(Z_c * Z * Z_c'));        
        
        %% IOC at-wall      
        C2_2 = flow.Ctaup(1:2,:);
        [Y,~,~] = care(flow.A',C2_2',flow.Bw*flow.Bw',Vwall);
        L = Y * C2_2' / Vwall;
        B2_2 = flow.Bbs(:,[1,7]);
        [X,~,~] = care(flow.A,B2_2,C1'*C1,Rp_wall);
        F = Rp_wall \ B2_2' * X;
        IOC_gam_taubbs(i,j) = trace(flow.Bw' * X * flow.Bw) + trace( inv(Rp_wall) * B2_2' * X*Y*X * B2_2);
        
        Z_a = [flow.A, -B2_2*F; L * C2_2, flow.A - B2_2 * F - L * C2_2];
        Z_b = [flow.Bw,zeros(size(flow.Bw,1),size(L*Vwall,2)) ; zeros(size(L*Vwall,1),size(flow.Bw,2)), L * sqrt(Vwall)]; % , zeros(size(flow.Bw,1),size(L*Vp_vel,2))];%; zeros(size(L*Vp05,1),size(flow.Bw,2)), L * Vp05];
        Z_c = [C1, zeros(size(C1,1),size(Rp_wall*F,2)); zeros(size(Rp_wall*F,1),size(C1,2)), -sqrt(Rp_wall) * F];
        Z = lyap(Z_a, Z_b * Z_b');
        IOC_rms_taubbs(i,j,:) = real(diag(Z_c * Z * Z_c'));


        disp(['i = ',num2str(i),', j = ',num2str(j),', toc = ',num2str(toc)])
    end
    toc
end
%%

% Uncontrolled flow
P_gam = sum(sum(P_gam_mat([1:end,2:end],[1:end,2:end]),2),1);
P_rms = squeeze(sum(sum(P_rms_output_mat([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1));

% Optimal Estimation / Actuation Everywhere in-flow
gam_u = sum(sum(OE_gam_u([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_u = squeeze(sum(sum(OE_rms_u([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

rms_u_forcingAE = squeeze(sum(sum(AE_rms_forcing([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1));% ./ P2_rms;

% Optimal Estimation / Actuation Everywhere at-wall
gam_tau = sum(sum(OE_gam_tau([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_tau = squeeze(sum(sum(OE_rms_tau([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

% Full information control everywhere in-flow
gam_z = sum(sum(FIC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_z = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

% Full information control everywhere at-wall
gam_bbs = sum(sum(FIC_gam_bbs([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_bbs = squeeze(sum(sum(FIC_rms_bbs([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

% Input output control in-flow
gam_u_z = sum(sum(IOC_gam_uz([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_u_z = squeeze(sum(sum(IOC_rms_uz([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

% Input output control at-wall
gam_tau_bbs = sum(sum(IOC_gam_taubbs([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_tau_bbs = squeeze(sum(sum(IOC_rms_taubbs([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

% Input output control mixed
gam_u_bbs = sum(sum(IOC_gam_ubbs([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
rms_u_bbs = squeeze(sum(sum(IOC_rms_ubbs([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;

%%
save('fig8','Re','Nx','Nboth','KX','KY',...
    'var_a','var_s','meas_loc_s','meas_loc_a', ...
    'P_gam','P_rms','gam_u','rms_u','gam_tau','rms_tau','gam_z','rms_z',...
    'gam_bbs','rms_bbs','gam_u_z','rms_u_z','gam_tau_bbs','rms_tau_bbs',...
    'gam_u_bbs','rms_u_bbs',...
    'P_gam_mat','P_rms_output_mat','FIC_rms_z','FIC_rms_bbs',...
    'FIC_gam_z','FIC_gam_bbs','IOC_gam_uz','IOC_rms_uz',...
    'IOC_gam_taubbs','IOC_rms_taubbs','OE_gam_u','OE_rms_u',...
    'AE_rms_forcing','OE_gam_tau','OE_rms_tau',...
    'IOC_gam_ubbs','IOC_rms_ubbs')


