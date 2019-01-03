function [sum_rate_opt_avg, norm_w_opt_avg, norm_v_opt_avg, EH_AP_opt_avg, EH_cyc_opt_avg] = ...
    fig2_convergence_function(P_S_dB, num_of_trial, r_U_th, r_D_th)

%% Add Solvers
% Extract SDPT3-4.0, SeDuMi_1_3, yalmip, etc. directly to the folder
pathCell = regexp(path, pathsep, 'split');
if (~any(strcmpi('.\yalmip', pathCell)))
    disp('Adding path...')
    addpath(genpath('./yalmip'));
    addpath(genpath('./SeDuMi_1_3'));
    addpath(genpath('./SDPT3-4.0'));
    addpath(genpath('./rician'));
end
clear pathCell

%% Initialization
% num_of_trial = 100;
max_iteration_per_trial = 30;

N_t = 3; % number of Tx antennas at BS
N_r = 3; % number of Rx antennas at BS

K = 0; % number of DLUs in inner zone
L = 1; % number of DLUs in outer zone

Q = 1; % number of ULUs in inner zone
P = 0; % number of ULUs in outer zone

% P_S_dB  = 20; % dowlink power (dBm)
P_S     = 10^(P_S_dB/10); % downlink power in miliwatt (mW)

% r_D_th      = 1; % data rate threshold of DL user (bps/Hz)
% r_U_th      = 1; % data rate threshold of UL user (bps/Hz)
threshold   = [r_D_th, r_U_th]*log(2); % convert to nats/sec/Hz

mu_S_dB     = -110; % degree of SI mitigation in dB (not dBm)
mu_U_dB     = 0; % no SI suppression
mu_S        = 10^(mu_S_dB/10);
mu_U        = 10^(mu_U_dB/10);

BW          = 10; % MHz
eta         = .1; % energy converter efficiency
Noise       = -80; % dBm
sigma_D     = sqrt(10^(Noise/10));
sigma_S     = sqrt(10^(Noise/10));
sigma       = [mu_S, mu_U, sigma_D, sigma_S, eta]; % 10^(Noise/10)*BW*10^6;

sumRateOptimalAll = zeros(num_of_trial, 1);
solutionOptimalAll = cell(num_of_trial, 1);
EH_AP_opt_all = zeros(num_of_trial, 1);
EH_cyc_opt_all = zeros(num_of_trial, 1);

%% Set users' positions
radius_cell = 20; % in meters (m)
radius_inner = 5;
radius_nearest_UE = 2;
std = 8; % standard deviation
ploss = 3; % path-loss exponent

Parameters = [radius_cell, radius_inner, radius_nearest_UE, std, ploss];

% Path-loss (but not used in this simulation)
% [D_H1k, D_H2l, D_G2p, D_G1q, D_f2_pq_kl, D_F1ql, ...
%     pos_DL_users, pos_UL_users] ...
%     = CreateLargeScaleFading(K, L, P, Q, Parameters);

%% Create Channels
channels_room = cell(num_of_trial, 9); % reserved rooms for channel realizations
v = sqrt(1); s = 1;
dist = {'Rice', v, s};
% The unit non-centrality parameter and the scale parameter
% K_factor = v^2/(2*s^2)

for ii = 1:1:num_of_trial
    % Path-loss setting used in the simulation
    D_H2l = eye(L)*sqrt(10^(-30/10)); % path-loss DL 30 dB
    D_G1q = eye(L)*sqrt(10^(-20/10)); % path-loss UL 20 dB (reciprocal)
    D_f2_pq_kl = ones(P+Q, K+L)*sqrt(10^(-20/10)); % path-loss CCI 20 dB
    
    channels_room{ii,1} ...
        = CreateSmallScaleFading(1, 1, N_t, L)*D_H2l; % hD
    channels_room{ii,2} ...
        = CreateSmallScaleFading(1, 1, N_t, Q)*D_G1q; % hU EH links
    channels_room{ii,3} ...
        = CreateSmallScaleFading(1, 1, N_r, Q)*D_G1q; % hU data links
    channels_room{ii,4} ...
        = CreateSmallScaleFading(1, 1, P+Q, K+L).*D_f2_pq_kl; % g_CCI
    channels_room{ii,5} ...
        = CreateSmallScaleFading(1, 1, N_t, N_r, dist); % G_SI
    channels_room{ii,6} ...
        = CreateSmallScaleFading(1, 1, 1, 1, dist); % g_SI at UL user
    channels_room{ii,7} ...
        = CreateSmallScaleFading(1, 1, N_r, L)*D_H2l; % DL half power
    channels_room{ii,8} ...
        = CreateSmallScaleFading(1, 1, N_t, L)*D_G1q; % UL half power
end

% Plot_Layout( radius_cell, radius_inner, pos_DL_users, pos_UL_users);
% h = waitbar(0, 'Please wait ...','Name','Percent of completion');

for ii = 1:1:num_of_trial
    % Assign channels
    hD          = channels_room{ii,1};
    hU_EH       = channels_room{ii,2};
    hU          = channels_room{ii,3};
    g_cci       = channels_room{ii,4};
    G_SI_BS     = channels_room{ii,5};
    g_SI_ULU    = channels_room{ii,6};
    hD_HalfP1   = channels_room{ii,7}; % un-used
    hU_HalfP1   = channels_room{ii,8}; % un-used
    % all channel relizations
    channels  = {hD, hU_EH, hU, g_cci, G_SI_BS, g_SI_ULU, hD_HalfP1, hU_HalfP1};
    
    % Notify completed percent
    %     percent = (ii-1)/num_of_running;
    %     waitbar(percent, h, [num2str(floor(percent*100)) ' % completed'])
    
    %% Proposed Algorithm
    [isBreak, sum_rate_opt, sum_rate_opt_chain, solution_opt, EH_AP_opt, EH_cyc_opt] ...
        = ProposedAlg(P_S, channels, threshold, sigma, max_iteration_per_trial, ii);
    
    %%
    %     disp('Convergence of the optimal sum-rate:')
    %     disp(sum_rate_opt_chain)
    
    sumRateOptimalAll(ii,1) = sum_rate_opt;
    solutionOptimalAll{ii,1} = solution_opt; % Note: this is a cell of (w,v,pU)
    %     EH_AP_opt_all(ii,1) = EH_AP_opt;
    %     EH_cyc_opt_all(ii,1) = EH_cyc_opt;
end

%% Print out the optimal average resutls
sum_rate_opt_avg = mean(sumRateOptimalAll);
EH_AP_opt_avg = mean(EH_AP_opt_all);
EH_cyc_opt_avg = mean(EH_cyc_opt_all);

norm_w_opt_sum  = zeros(1,1);
norm_v_opt_sum  = zeros(1,1);

for jj = 1:num_of_trial
    norm_w_opt_sum = norm_w_opt_sum + norm(solutionOptimalAll{jj,1}{1,1})^2;
    norm_v_opt_sum = norm_v_opt_sum + norm(solutionOptimalAll{jj,1}{1,2})^2;
end

norm_w_opt_avg = norm_w_opt_sum/num_of_trial;
norm_v_opt_avg = norm_v_opt_sum/num_of_trial;
% close(h)
end