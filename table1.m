% This file generates displays the data for Figure 10 for
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

clear 

load('fig5/fig5.mat')

% select kx up to 0.5
Kximax = 3;
% select ky up to 6
Kyjmax = 10;

% calculate H2-norm^2
P_gam = sum(sum(P_gam_mat([1:Kximax,2:Kximax],[1:Kyjmax,2:Kyjmax]),2),1);
gam2_u = (sum(sum(OE_gam_u([1:Kximax,2:Kximax],[1:Kyjmax,2:Kyjmax]),2),1) / P_gam);
gam2_z = (sum(sum(FIC_gam_z([1:Kximax,2:Kximax],[1:Kyjmax,2:Kyjmax]),2),1) / P_gam);
gam2_uz = (sum(sum(IOC_gam_z([1:Kximax,2:Kximax],[1:Kyjmax,2:Kyjmax]),2),1) / P_gam);

disp('AE:')
disp(real(gam2_u))
disp('ME:')
disp(real(gam2_z))
disp('IO:')
disp(real(gam2_uz))

