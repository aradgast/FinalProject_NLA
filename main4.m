clc
clear
close all

% Parameters
lamda = 1;
W = 64 * lamda;
delta = lamda / 10;
N = W / delta;
angle = pi / 2;
alpha = 1;

%%

A = create_steering_mat(lamda, W, angle, alpha);
[U,S,V] = svd(A);
t_orig = original_transmitter_weights(lamda, W, angle, alpha);
error_delta_range = logspace(-15, 0, 50);

% solve with full SVD
[~, error_list_full_SVD] = find_t_and_error(A, U, S, V, t_orig, error_delta_range);

condition_num_A = cond(A);
disp('kappa value:')
disp(num2str(condition_num_A,5))

% solve with truncated SVD
tau_list = [10^-2, 10^-5, 10^-8];
error_list_truncated_SVD = zeros(length(tau_list), length(error_delta_range));
for tau_idx=1:length(tau_list)
    tau = tau_list(tau_idx);
    rank = sum(diag(S) > tau);
    U_t = U(:, 1:rank);
    V_t = V(:, 1:rank);
    S_t = S(1:rank, 1:rank);
    [~, error_list_tmp] = find_t_and_error(A, U_t, S_t, V_t, t_orig, error_delta_range);
    error_list_truncated_SVD(tau_idx, :) = error_list_tmp;
end

%%
figure(1)
plot(error_delta_range, error_list_full_SVD,'LineWidth',2)
hold on
plot(error_delta_range, error_list_truncated_SVD,'LineWidth',2)
legend("full SVD", "Trancated, \tau = 10^{-2}", "Trancated, \tau = 10^{-5}", "Trancated, \tau = 10^{-8}")
title("4.m - Error as a function of \delta")
xlabel("\delta")
ylabel("MSE")
set(gca,'xscale','log','YScale','log')
grid on

hold off

condition_num_A = cond(A);
disp('kappa value:')
disp(num2str(condition_num_A,5))

%% n - solve with Different Algo.1

tau_list = [10^-2, 10^-5, 10^-8];
gamma = 10;
error_list_algosolve_SVD = zeros(length(tau_list), length(error_delta_range));
algo_num = 2; %1 or 2
monte_carlo_num = 10; %change the mc number to 1 inside find_t_and_error_func

error_list_algosolve_SVD = zeros(2, length(tau_list), length(error_delta_range));
for t=1:size(tau_list,2)
    tau = tau_list(1,t);
    if algo_num == 1
        tau_r_or_B_0 = 0.1; %for tau_r
    elseif algo_num == 2
        tau_r_or_B_0 = 1; %for B_0
    else 
        disp("Wrong Algo_num values")
    end
    for mc=1:monte_carlo_num
        [U_LR_hat,S_LR_hat,V_LR_hat, Rank_l] = Improve_accuracy_for_LS_solve(algo_num, A, gamma, tau, tau_r_or_B_0);
        [~ ,error_list_tmp] = find_t_and_error(A,U_LR_hat,S_LR_hat,V_LR_hat,t_orig, error_delta_range);
        error_list_algosolve_SVD(algo_num, t,:) = error_list_algosolve_SVD(algo_num, t, :) + reshape(error_list_tmp, [], 1, length(error_list_tmp)) / monte_carlo_num;
    end
end

%plot error found with Algo num
figure(2)
plot(error_delta_range, error_list_full_SVD,'LineWidth',2)
hold on
plot(error_delta_range, squeeze(error_list_algosolve_SVD(algo_num, :, :)),'LineWidth',2)
title("4.n - MSE as a function of \delta")
if algo_num == 1
    subtitle("Using ILRA Alg., with \gamma = " + num2str(gamma))
else if algo_num == 2
    subtitle("Using UILRA Alg., with B_{0} = " + num2str(tau_r_or_B_0) + " and \gamma = " + num2str(gamma))
end
end
xlabel("\delta")
ylabel("MSE")
set(gca,'xscale','log','YScale','log')
grid on
legend("full SVD", '\tau = 10^-2','\tau = 10^-5','\tau = 10^-8')
hold off



%% PINV solution as a reference
monte_carlo_num = 100;
r = A * t_orig; %Nx1
N = size(A, 1);
A_pinv = pinv(A);

error_list_pinv = zeros(1, length(error_delta_range));
t_tild_list = zeros(size(r,1), length(error_delta_range));

for d = 1:length(error_delta_range)
    error_delta = error_delta_range(1, d);
    for mc=1:monte_carlo_num
        teta = unifrnd(0, 2 * pi, N, 1);
        r_plagued = r + error_delta * (abs(r).*exp(1i * teta));
        t_plagued = A_pinv * r_plagued; %  N x 1
        t_tild_list(:,d) = t_plagued;
        error_list_pinv(1,d) = error_list_pinv(1,d) + norm(t_orig-t_plagued,2) / (length(t_orig) * monte_carlo_num);

    end
end
%%
%compare all:L fill SVD, truncted SVD, algo1 SVD, algo2 SVD
figure(3)
plot(error_delta_range, error_list_full_SVD,'LineWidth', 2,'Color',[0,0,0])
hold on
plot(error_delta_range, error_list_pinv, 'LineWidth', 2) % Access marker using curly braces
SVDfor_plot = [error_list_truncated_SVD; squeeze(error_list_algosolve_SVD(1, :, :)); squeeze(error_list_algosolve_SVD(2, :, :))];
marker = {'-*', '-x', '-|'};
idx_marker = 0;
for idx = 1:size(SVDfor_plot, 1)
    if mod(idx, 3) == 1
        idx_marker = idx_marker + 1;
    end
    plot(error_delta_range, SVDfor_plot(idx, :), marker{idx_marker}, 'LineWidth', 2)

end

legend("full SVD", "PINV", ...
    'Trancated, \tau = 10^-2', 'Trancated, \tau = 10^-5', 'Trancated, \tau = 10^-8', ...
    'ILRA \tau = 10^-2', 'ILRA \tau = 10^-5', 'ILRA \tau = 10^-8', ...
    'UILRA \tau = 10^-2', 'UILRA \tau = 10^-5', 'UILRA \tau = 10^-8')
title("4.n - Error as a function of \delta (Comparing all)")

xlabel("\delta")
ylabel("MSE")
set(gca,'xscale','log','YScale','log')
grid on
hold off
    
    



