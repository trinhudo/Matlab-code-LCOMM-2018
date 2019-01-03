close all
clear all
clc

%% Parameters
num_of_trial = 50; % number of trials
P_S_dB = 10:5:35; % dowlink power budget (dBm)
r_U_th = 1; % bps/Hz
r_D_th = 1; % bps/Hz

sum_rate_opt_avg = zeros(1,length(P_S_dB));
norm_w_opt_avg   = zeros(1,length(P_S_dB));
norm_v_opt_avg   = zeros(1,length(P_S_dB));
EH_AP_opt_avg    = zeros(1,length(P_S_dB));
EH_cyc_opt_avg   = zeros(1,length(P_S_dB));

%% Get optimal results
for ii = 1:length(P_S_dB)
    [sum_rate_opt_avg(ii), norm_w_opt_avg(ii), norm_v_opt_avg(ii), ...
        EH_AP_opt_avg(ii), EH_cyc_opt_avg(ii)] = ...
    fig2_convergence_function(P_S_dB(ii), num_of_trial, r_U_th, r_D_th);
end

%% Print optimal results
format short
disp('--- Final Optimal Results ---')

disp('The transmit power budget (mW) is:')
disp(10.^(P_S_dB./10))

disp('The optimal average sum-rate (bps/Hz) is:')
disp(sum_rate_opt_avg)

disp('The optimal average ||w||^2 is:')
disp(norm_w_opt_avg)

disp('The optimal average ||v||^2 is:')
disp(norm_v_opt_avg)

disp('The harvested energy (mW) from H-AP at UL user is:')
disp(EH_AP_opt_avg)

disp('The energy recycled (mW) at UL user is:')
disp(EH_cyc_opt_avg)