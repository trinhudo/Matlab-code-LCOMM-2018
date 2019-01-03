function [isBreak, sum_rate_opt, sum_rate_opt_chain, solution_opt, EH_AP_opt, EH_cyc_opt] ...
    = ProposedAlg(P_S,  channels, thresholds, sigma, max_iteration_per_trial, ii_num_of_trial)
% Inputs
% P_S: power budget at Source in mW
% channels: all channel realizations 
% thresholds: UL and DL data rate threshold in nats/sec/Hz
% sigma: all noises and eta, namely: mu_S = sigma(1), mu_U = sigma(2),
% sigma_D = sigma(3), sigma_S = sigma(4), eta = sigma(5)

% Outputs
% solution_opt: the optimal solution of (w,v,pU)

L = size(channels{1,1}, 2);
Q = size(channels{1,2}, 2);
[N_t, ~] = size(channels{1,5});

%% Loop for optimization
GetBreak = 0;
if (GetBreak)
    isBreak = 1;
    global_OptValue = 0;
    global_OptValueChain = [];
    global_OptValueChain_sum = 0;
    global_UplinkRate_PerGroupPerUser = 0;
    global_DownlinkRate_PerGroupPerUser = 0;
else
    isBreak = 0;
    disp('Run iterative algorithm...');
    nn = 0;
    OptimalValue_current = 0;
    OptValueChain = [];
    RunningMode = 'Init';
    preStatus = {'', 0, RunningMode};
    InitPoint = -1;
    
    w_current = 1*(randn(N_t,L)+1i*randn(N_t,L));
    v_current = 1*(randn(N_t,1)+1i*randn(N_t,1));
%     v_current = zeros(N_t,1); % for the case of no energy beam
    pU_current = 1*rand(1, Q);
    X_current = {w_current, v_current, pU_current};
    z_current = (norm(channels{1,1}'*v_current)^2 ...
        + pU_current^2*abs(channels{1,4})^2 ...
        + sigma(3)^2)/((real(channels{1,1}'*w_current))^2);
    
    while (nn < max_iteration_per_trial)
        fprintf('P_S_dB = %d dB, %d-th trial, %d-th iteration \n', ...
            10*log10(P_S), ii_num_of_trial, nn+1);
        % Get optimal solution in the current iteration
        [OptimalValue, X_next, EH_AP_opt, EH_cyc_opt, z_next, Status] ...
            = Get_optSolutionPerIteration(P_S, channels, ...
            thresholds, sigma, X_current, z_current, preStatus);
        % Update the optimal solution for the next iteration
        if (isFeasible(OptimalValue_current, OptimalValue, ...
                InitPoint, preStatus, Status, RunningMode))
            X_current = X_next;
            z_current = z_next;
            ww = trace(X_current{1,1}'*X_current{1,1});
            vv = trace(X_current{1,2}'*X_current{1,2});
            SumPbs = ww + vv; % sum of weights <= power budget
            pU_sqr = X_current{1,3}.^2; % transmit power of UL user
            %
            if ((~isempty(findstr(Status,'Successfully'))) ...
                    && OptimalValue > 0)
                RunningMode = '';
            end
            %
            if (~isempty(findstr(RunningMode,'Init')))
                InitPoint = InitPoint+1;
            end
            preStatus = {Status, OptimalValue, RunningMode};
            % Check convergence
            if (checkConvergence(OptValueChain, OptimalValue) ...
                    && (OptimalValue>0))
                break;
            else
                if ((isempty(findstr(RunningMode,'Init'))) ...
                        && OptimalValue > 0)
                    OptValueChain = [OptValueChain OptimalValue];
                else
                    OptValueChain = [];
                end
                OptimalValue_current = OptimalValue;
            end
        else
            disp('--> keeping the latest feasible point');
            OptValueChain = [OptValueChain OptimalValue_current];
            break;
        end
        nn = nn + 1;
    end
    
    solution_opt = X_current; % (w,v,pU)
    sum_rate_opt = OptimalValue_current/log(2); % convert to bps/Hz
    sum_rate_opt_chain = OptValueChain/log(2); % convert to bps/Hz
    
    fprintf('The convergence of the %d-th trial is: \n', ii_num_of_trial)
    disp(sum_rate_opt_chain)
    
    fprintf('The optimal solution of the %d-th trial is: \n', ii_num_of_trial)
    disp(solution_opt)
end
end
%%
function [check] = isFeasible(OptimalValue_current, OptimalValue, ...
    InitPoint, preStatus, Status, RunningMode)
% check if the current point is feasible and the next point is infeasible
check = 1;
if ((isempty(findstr(Status,'Successfully'))) ...
        && OptimalValue_current>0 ...
        && (~isempty(findstr(preStatus{1,1},'Successfully'))) ...
        && (isempty(findstr(RunningMode,'Init'))))
    disp('Good job!');
    check = 0;
end
%
check_time = 1;
check_time_next = 1;
check3 = 1;
if ((OptimalValue_current > OptimalValue) ...
        && (InitPoint > 0) && (OptimalValue > 0))
    check3 = 0;
end
check = (check_time && check_time_next && check3 && check);
end
%%
function [check] = checkConvergence(OptValueChain, OptimalValue)
check = 0;
if (length(OptValueChain)>10)
    if (abs(OptimalValue-OptValueChain(length(OptValueChain))) < 10^-2)
        check = 1;
    end
    if (abs(OptimalValue-OptValueChain(length(OptValueChain)-5)) < 0.01)
        check = 1;
    end
end
end