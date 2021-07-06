% This file generates the support files for Figure 11 of
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

clear

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

% set alpha
alphas = [1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
labels = {'n5','n4','n3','n2','n1','1','p1','p2'};

alphas = [1e-1,1e0,1e1,1e2];
labels = { 'n1','1','p1','p2'};

k = 0;

for alpha = alphas
    
    k = k + 1;
    
    Rp_vel = (eye(3) * alpha)^2;
    Vp_vel = (eye(3) * 1e-4)^2;
    
    
    %% simulate linear model at each wavenumber pair
    
    
    P_rms_output_mat = zeros(length(KX),length(KY),  3*(flow.N_out + 1));
    P_gam_mat = zeros(length(KX),length(KY));
    
    FIC_gam_z = zeros(length(KX),length(KY));
    FIC_rms_z = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + 3);
    
    IOC_gam_z = zeros(length(KX),length(KY));
    IOC_rms_z = zeros(length(KX),length(KY),  3*(flow.N_out + 1) + 3);
    
    %%
    var_a = 0.02;
    var_s = 0.02;
    
    meas_loc_s = 0.68;%1 - 0.259; % incorrect
    meas_loc_a = 0.68;%1 - 0.325; % incorrect
    Bint_dz = exp( -((flow.geom.x - meas_loc_a) / var_a ) .^2 );
    Bint_dz = [Bint_dz, zeros(size(Bint_dz)), zeros(size(Bint_dz)); ...
        zeros(size(Bint_dz)), Bint_dz, zeros(size(Bint_dz)); ...
        zeros(size(Bint_dz)), zeros(size(Bint_dz)), Bint_dz];
    Cint_u = ( diag(flow.w_out) * exp( -((flow.z_out - meas_loc_s) / var_s ) .^2 ) ).';
    Cint_u = [Cint_u, zeros(size(Cint_u)), zeros(size(Cint_u)); ...
        zeros(size(Cint_u)), Cint_u, zeros(size(Cint_u)); ...
        zeros(size(Cint_u)), zeros(size(Cint_u)), Cint_u];
    
    
    %%
    for i = 1:length(KX)
        tic
        kx = KX(i);
        parfor j = 1:length(KY)
            tic
            % linear model at this wavenumber pair
            flow = channelOSS.StateSpace(kx,KY(j),Re,Nx,'turbulent','top',Nboth);

            C1 = flow.Cw;
            
            %% Uncontrolled
            Z2 = lyap(flow.A, flow.Bw * flow.Bw');
            P_gam_mat(i,j) = real(trace(C1 * Z2 * C1'));
            P_rms_output_mat(i,j,:) = real(diag(C1 * Z2 * C1'));
            %% AE
            C2 = Cint_u * flow.C;
            [Y,~,~] = care(flow.A',C2',flow.Bw*flow.Bw',Vp_vel);
            %% ME
            B2_1 = flow.B * Bint_dz;
            [X,~,~] = care(flow.A,B2_1,C1'*C1,Rp_vel);
            FIC_gam_z(i,j) = real(trace(flow.Bw' * X * flow.Bw));
            
            F = Rp_vel \ B2_1' * X;
            Z_a = flow.A - B2_1 * F;
            Z_b = flow.Bw;
            Z_c = [C1;-sqrt(Rp_vel) * F];
            X2 = lyap(Z_a, Z_b * Z_b');
            FIC_rms_z(i,j,:) = real(diag(Z_c * X2 * Z_c'));
            %% IOC
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
    %%
    P_gam = sum(sum(P_gam_mat([1:end,2:end],[1:end,2:end]),2),1);
    P_rms = squeeze(sum(sum(P_rms_output_mat([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1));
    
    %%
    gam_z = sum(sum(FIC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
    rms_z = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;
    
    gam_uz = sum(sum(IOC_gam_z([1:end,2:end],[1:end,2:end]),2),1) / P_gam;
    rms_uz = squeeze(sum(sum(IOC_rms_z([1:end,2:end],[1:end,2:end],1:Nboth+1),2),1)) ./ P_rms;
    
    
    save(['alpha_',labels{k}],'P_gam','P_rms','gam_z','rms_z','gam_uz','rms_uz','Rp_vel')
    
end

toc
