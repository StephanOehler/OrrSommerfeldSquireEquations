% This file generates the support files for Figure 10 of
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

% measurement and estimation locations

Rp_vel = (1e-4)^2;
Vp_vel = (1e-4)^2;


%% preallocate memory

P_rms_output_mat = zeros(length(KX),length(KY),  3*(flow.N_out + 1));
P_gam_mat = zeros(length(KX),length(KY));

FIC_gam_z = zeros(length(KX),length(KY));
FIC_rms_z = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + 1);


%%
% sensor and actuator gaussian width
var_a = 0.02;
var_s = 0.02;

% location of sensor plane

meas_loc_s = 0.68;
meas_loc_a = 0.68;



Cint_u_simple = ( diag(flow.w_out) * exp( -((flow.z_out - meas_loc_s) / var_s ) .^2 ) ).';
Cint_u{1} = [Cint_u_simple, zeros(size(Cint_u_simple)), zeros(size(Cint_u_simple))];
Cint_u{2} = [zeros(size(Cint_u_simple)), Cint_u_simple, zeros(size(Cint_u_simple))];
Cint_u{3} = [zeros(size(Cint_u_simple)), zeros(size(Cint_u_simple)), Cint_u_simple];

Bint_dz_simple = exp( -((flow.geom.x - meas_loc_a) / var_a ) .^2 );
Bint_dz{1} = [Bint_dz_simple; ...
    zeros(size(Bint_dz_simple)); ...
    zeros(size(Bint_dz_simple))];
Bint_dz{2} = [zeros(size(Bint_dz_simple)); ...
    Bint_dz_simple; ...
    zeros(size(Bint_dz_simple))];
Bint_dz{3} = [zeros(size(Bint_dz_simple)); ...
    zeros(size(Bint_dz_simple)); ...
    Bint_dz_simple];


%%

names = {'xu','yv','zw'};


for xyz = 1:3
    
    for i = 1:length(KX)
        tic
        kx = KX(i);
        for j = 1:length(KY)
            tic
            %% linear model at this wavenumber pair
            flow = channelOSS.StateSpace(kx,KY(j),Re,Nx,'turbulent','top',Nboth);
            C1 = flow.Cw;
            %% Uncontrolled flow
            Z2 = lyap(flow.A, flow.Bw * flow.Bw');
            P_gam_mat(i,j) = real(trace(C1 * Z2 * C1'));
            P_rms_output_mat(i,j,:) = real(diag(C1 * Z2 * C1'));
            %% ME / FIC
            B2_1 = flow.B * Bint_dz{xyz};
            [X,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel);
            FIC_gam_z(i,j) = real(trace(flow.Bw' * X * flow.Bw));
            
            F = Rp_vel \ B2_1' * X;
            Z_a = flow.A - B2_1 * F;
            Z_b = flow.Bw;
            Z_c = [C1;-sqrt(Rp_vel) * F];
            X2 = lyap(Z_a, Z_b * Z_b');
            FIC_rms_z(i,j,:) = real(diag(Z_c * X2 * Z_c'));        
            %%
            disp(['i = ',num2str(i),', j = ',num2str(j),', toc = ',num2str(toc)])
        end
        toc
    end
    %%
    P_gam = sum(sum(P_gam_mat([1:end,2:end],[1:end,2:end]),2),1);
    P_rms = squeeze(sum(sum(P_rms_output_mat([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1));
 
    gam_z = sum(sum(FIC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
    rms_z = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;
    
    save(['components_',names{xyz}])    
    
end
