% This file generates the support files for Figure 10 of
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

% This code uses Input Output control to mimic Actuation everywhere but
% only with certain components active (either x, y or z).

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


% cases
names = {'x','y','z'}; 
 
for xyz = 1:length(names)

    
    %% actuator cost and sensor noise 
    Rp_vel = (eye(Nx+1) * 1e-4)^2;
    Vp_vel = (eye(3) * 1e-4)^2;
    
    
    %% preallocate memory
    P_rms_output_mat = zeros(length(KX),length(KY),  3*(flow.N_out + 1));
    P_gam_mat = zeros(length(KX),length(KY));
    
    IOC_gam_z = zeros(length(KX),length(KY));
    IOC_rms_z = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + Nx+1);
    
    %%
    % sensor and actuator gaussian width
    var_a = 0.02;
    var_s = 0.02;
    
    % location of sensor plane
    meas_loc_s = 0.68;% 1 - 0.32
    
    % create sesnor matrix
    Cint_u = ( diag(flow.w_out) * exp( -((flow.z_out - meas_loc_s) / var_s ) .^2 ) ).';
    Cint_u = [Cint_u, zeros(size(Cint_u)), zeros(size(Cint_u)); ...
        zeros(size(Cint_u)), Cint_u, zeros(size(Cint_u)); ...
        zeros(size(Cint_u)), zeros(size(Cint_u)), Cint_u];
    
    
     % control everywhere matrix
     Bint_dz{1} = [eye(Nx+1); ...
     zeros(Nx+1); ...
     zeros(Nx+1)];
     Bint_dz{2} = [zeros(Nx+1); ...
     eye(Nx+1); ...
     zeros(Nx+1)];
     Bint_dz{3} = [zeros(Nx+1); ...
     zeros(Nx+1); ...
     eye(Nx+1)]; 

    %%
    for i = 1:length(KX)
        tic
        kx = KX(i);
        parfor j = 1:length(KY)
            tic
            %% linear model at this wavenumber pair
            flow = channelOSS.StateSpace(kx,KY(j),Re,Nx,'turbulent','top',Nboth);
            C1 = flow.Cw;
            
            %% Calculate Energy of uncontrolled flow
            Z2 = lyap(flow.A, flow.Bw * flow.Bw');
            P_gam_mat(i,j) = real(trace(C1 * Z2 * C1'));
            P_rms_output_mat(i,j,:) = real(diag(C1 * Z2 * C1'));
            %% calculate energy of controlled case
            C2 = Cint_u * flow.C;
            [Y,~,~] = care(flow.A',C2',flow.Bw*flow.Bw',Vp_vel);
            
            B2_1 = flow.B * Bint_dz{xyz};
            [X,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel);
            
            F = Rp_vel \ B2_1' * X;
            L = Y * C2' / Vp_vel;
            
            IOC_gam_z(i,j) = trace(flow.Bw' * X * flow.Bw) + trace( inv(Rp_vel) * B2_1' * X*Y*X * B2_1);
            
            Z_a = [flow.A, -B2_1*F; L * C2, flow.A - B2_1 * F - L * C2];
            Z_b = [flow.Bw,zeros(size(flow.Bw,1),size(L*Vp_vel,2)) ; zeros(size(L*Vp_vel,1),size(flow.Bw,2)), L * sqrt(Vp_vel)]; % , zeros(size(flow.Bw,1),size(L*Vp_vel,2))];%; zeros(size(L*Vp05,1),size(flow.Bw,2)), L * Vp05];
            Z_c = [C1, zeros(size(C1,1),size(Rp_vel*F,2)); zeros(size(Rp_vel*F,1),size(C1,2)), -sqrt(Rp_vel) * F];
            Z = lyap(Z_a, Z_b * Z_b');
            IOC_rms_z(i,j,:) = real(diag(Z_c * Z * Z_c'));
            
            disp(['i = ',num2str(i),', j = ',num2str(j),', toc = ',num2str(toc)])
        end
        toc
    end
    %% calculate H2 norms.
    P_gam = sum(sum(P_gam_mat([1:end,2:end],[1:end,2:end]),2),1);
    gam_z = sum(sum(IOC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
    
    % save data
    save(['components_',names{xyz}])


end
