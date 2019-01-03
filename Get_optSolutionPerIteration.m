function [sum_rate_opt, solution_opt, EH_AP_opt, EH_cyc_opt, z, Status] = ...
    Get_optSolutionPerIteration(P_S, channels, thresholds, sigma, ...
    X_current, z_current, preStatus)
% Inputs
% P_S: power budget at Source in mW
% channels: all channel realizations 
% thresholds: UL and DL data rate threshold in nats/sec/Hz
% sigma: noise spectral in dBm (not dBm/Hz)

% Outputs
% sum_rate_opt: optimal sum rate (nats/sec/Hz)
% solution_opt = (w_opt, v_opt, pU_opt)
% EH_AP: energy harvested from H-AP
% EH_cyc: energy harvested from recycling

if (~isempty(findstr(preStatus{1,3},'Init')))
    disp('Initlizing... ');
end

%% Parameters
L           = size(channels{1,1}, 2);
Q           = size(channels{1,4}, 2);
[N_t, ~]    = size(channels{1,5});

% Assign channels, SI, and AWGN
hD          = channels{1,1};
hU_EH       = channels{1,2};
hU          = channels{1,3};
g_cci       = channels{1,4};
G_SI_BS     = channels{1,5};
g_SI_UL     = channels{1,6};

mu_S        = sigma(1);
mu_U        = sigma(2);
sigma_D     = sigma(3);
sigma_S     = sigma(4);
eta         = sigma(5);

r_D_th      = thresholds(1); % nats/sec/Hz
r_U_th      = thresholds(2);

% Extract current solution
w_current   = X_current{1,1};
v_current   = X_current{1,2};
pU_current  = X_current{1,3};

% Initialize variables
w_next      = sdpvar(N_t, L, 'full','complex');
v_next      = sdpvar(N_t, 1, 'full','complex');
% v_next      = zeros(N_t,1); % for the case of no energy beam
pU_next     = sdpvar(1, Q, 'full','real');
varphi      = sdpvar(1, 1, 'full','real');

% Add objective function and constraints
% obj = 0;
cons = [];

%% Constraints
for iQ = 1:1:Q
    cons = [cons, pU_next(iQ) >= 10^(-7)]; 
    % additional contraint to make sure positve harvested energy
end
cons    = [cons, cone([w_next; v_next], sqrt(P_S))];

z_next  = sdpvar(1, 1, 'full','real');
R_SD    = sdpvar(1, 1, 'full','real'); % Source-DL user data rate
Phi     = sdpvar(1, 1, 'full','real');

Psi     = log(1+1/z_current) + 1/(z_current+1);
Xi      = 1/(z_current*(z_current+1));
R_SD    = Psi - Xi*z_next;
Phi     = real(hD'*w_current)*(2*real(hD'*w_next)-real(hD'*w_current));

cons    = [cons, Phi>=0];

cons    = [cons, cone([hD'*v_next, pU_next*g_cci, ...
    sigma_D, 0.5*(z_next-Phi)], 0.5*(z_next+Phi))];

cons    = [cons, z_next >= 10^(-7)];

cons    = [cons, real(R_SD)>=0];

if (isempty(findstr(preStatus{1,3},'Init')))
    cons = [cons, real(R_SD) - r_D_th >= 0]; % for initial point
else
    cons = [cons, real(R_SD)-r_D_th>=varphi]; % for the convex problem
end

R_US    = sdpvar(1, 1, 'full','real'); % UL user-Source data rate
t_next  = sdpvar(1, 1, 'full','real');

Theta   = mu_S*G_SI_BS'*w_current*w_current'*G_SI_BS ...
    + mu_S*G_SI_BS'*v_current*v_current'*G_SI_BS ...
    + sigma_S^2*eye(size(G_SI_BS,2));

Theta_d  = inv(Theta) - inv(Theta + pU_current^2*hU*hU');

gamma_US = real(pU_current^2*hU'*inv(Theta)*hU);

cons = [cons, cone([sqrt(mu_S)*w_next'*G_SI_BS*Theta_d^(0.5), ...
    sqrt(mu_S)*v_next'*G_SI_BS*Theta_d^(0.5), ...
    sigma_S*sqrt(trace(Theta_d)), pU_next*hU'*Theta_d^(0.5), ...
    0.5*(t_next-1)], 0.5*(t_next+1))];

cons = [cons, t_next >= 10^(-7)];

R_US = log(1+gamma_US)-gamma_US ...
    + 2*pU_current*hU'*inv(Theta)*hU*pU_next - t_next;

cons = [cons, real(R_US)>=0];

if (isempty(findstr(preStatus{1,3},'Init')))
    cons = [cons, real(R_US) - r_U_th >= 0 ];
else
    cons = [cons, real(R_US)-r_U_th>=varphi];
end

pU_EH = sdpvar(1, 1, 'full','real'); % total harvested energy at UL user

pU_EH = eta*( 2*real(w_current'*hU_EH*hU_EH'*w_next) ...
    - norm(hU_EH'*w_current)^2 ...
    + 2*real(v_current'*hU_EH*hU_EH'*v_next) ...
    - norm(hU_EH'*v_current)^2 ...
    + mu_U*norm(g_SI_UL)^2*(2*pU_current*pU_next-pU_current^2));

cons = [cons, real(pU_EH) >= 10^(-7)];

cons = [cons, cone([pU_next, 0.5*(pU_EH-1)], 0.5*(pU_EH+1))];

%% Call solver
if (isempty(findstr(preStatus{1,3},'Init')))
    cons    = [cons, (real(R_SD) + real(R_US))>=varphi]; 
    obj     = varphi; % for initial point
else
    cons    = [cons, varphi>=0]; 
    obj     = varphi; % for the convex problem
end

myops       = sdpsettings('solver','sdpt3','verbose',0); % using SDPT-3
diagnotics  = solvesdp(cons, -obj, myops);

%% Get optimal results
sum_rate_opt = double(obj);
solution_opt = {double(w_next), double(v_next), double(pU_next)}; % cells
z            = double(z_next);
pU_EH_opt    = double(pU_EH);
DL_rate_opt  = double(real(R_SD));
UL_rate_opt  = double(real(R_US));

EH_AP_opt    = eta*(norm(hU_EH'*w_next)^2 + norm(hU_EH'*v_next)^2);
EH_cyc_opt   = eta*(mu_U*(pU_next^2)*(norm(g_SI_UL)^2));

Status       = diagnotics.info;
end