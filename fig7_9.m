% This file generates displays the data for figures 7 and 9 in
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735


% clear workspace
clear; close all

%
load('fig8/fig8.mat')

% select indicies
idx_iii = [1:Nboth+1;Nboth+2:2*Nboth+2;2*Nboth+3:3*Nboth+3];


%% Uncontrolled flow
data_control = squeeze(sum(sum(P_rms_output_mat([1:end,2:end],[1:end,2:end],:),2),1));
Pall = sum(data_control);
Puvw = sum(data_control(idx_iii),2);

disp('Uncontrolled flow')
disp(Puvw / Pall * 100)

%% AE / FAC
data_control = squeeze(sum(sum(OE_rms_u([1:end,2:end],[1:end,2:end],:),2),1)) ;
FACall = real(sum(data_control));
FACuvw = real(sum(data_control(idx_iii),2));

data_control_FAC = squeeze(sum(sum(AE_rms_forcing([1:end,2:end],[1:end,2:end],:),2),1)) ;
FACqforce = real(sum(data_control_FAC(idx_iii),2));

% select the desired data
disp('Actuation everywhere')
disp(real([FACall;FACuvw]) / Pall * 100)

% display the desired data
disp('Actuation everywhere - distribution of forcing')
disp(FACqforce / sum(FACqforce) * 100)

%% FIC / ME
% select the desired data
data_control = squeeze(sum(sum(FIC_rms_z([1:end,2:end],[1:end,2:end],:),2),1)) ;
FICall = real(sum(data_control));
FICuvw = real(sum(data_control(idx_iii),2));

% display the desired data
disp('Measurements Everywhere Control')
disp(real([FICall;FICuvw]) / Pall * 100)
disp('Measurements Everywhere Control - distribution of forcing')
disp(data_control(end-2:end) ./ sum(data_control(end-2:end)) *100)

%% IOC
% select the desired data
data_control = squeeze(sum(sum(IOC_rms_uz([1:end,2:end],[1:end,2:end],:),2),1)) ;
IOCall = real(sum(data_control));
IOCuvw = real(sum(data_control(idx_iii),2));

% display the desired data
disp('Input Output Control')
disp(real([IOCall;IOCuvw]) / Pall * 100)
disp('Input Output Control - distribution of forcing')
disp(data_control(end-2:end) ./ sum(data_control(end-2:end)) *100)






