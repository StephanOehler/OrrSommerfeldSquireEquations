% This file generates displays the data for Figure 10 for
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

names = {'x','y','z'}; 

IOC_vec_uvw = zeros(1,3);
IOC_vec_u = IOC_vec_uvw;
IOC_vec_v = IOC_vec_uvw;
IOC_vec_w = IOC_vec_uvw;

N_out = 200;
for xyz = 1:3

load(['fig10_11/components_',names{xyz}],'gam_z','IOC_rms_z','P_gam')

rms_temp = squeeze(sum(sum(IOC_rms_z([1:end,2:end],[1:end,2:end],:))));
IOC_vec_uvw(xyz) = sum(rms_temp(1:3*N_out+3))/P_gam;
IOC_vec_u(xyz) = sum(rms_temp(1:N_out+1))/P_gam;
IOC_vec_v(xyz) = sum(rms_temp(N_out+2:2*N_out+2))/P_gam;
IOC_vec_w(xyz) = sum(rms_temp(2*N_out+3:3*N_out+3))/P_gam;

end

disp('Fig 10: AE')
disp([real(IOC_vec_uvw);real(IOC_vec_u);real(IOC_vec_v);real(IOC_vec_w)]);

%%

names = {'xu','yv','zw'}; 

ME_vec_uvw = size(1,3);
ME_vec_u = ME_vec_uvw;
ME_vec_v = ME_vec_uvw;
ME_vec_w = ME_vec_uvw;

N_out = 200;
for xyz = 1:3 

load(['fig10_11/components_',names{xyz}],'gam_z','FIC_rms_z','P_gam')

rms_temp = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],1:(3*N_out+3)))));

ME_vec_uvw(xyz) = sum(rms_temp(1:3*N_out+3))/P_gam;
ME_vec_u(xyz) = sum(rms_temp(1:N_out+1))/P_gam;
ME_vec_v(xyz) = sum(rms_temp(N_out+2:2*N_out+2))/P_gam;
ME_vec_w(xyz) = sum(rms_temp(2*N_out+3:3*N_out+3))/P_gam;
end

disp('Fig 10: ME')
disp([real(ME_vec_uvw);real(ME_vec_u);real(ME_vec_v);real(ME_vec_w)])

