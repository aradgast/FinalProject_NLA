clc
clear
close all

% Parameters
lamda = 1;
W = 4 * lamda;
delta = lamda / 10;
N = W / delta;
angle = pi / 2;
alpha = 1;

A = create_steering_mat(lamda, W, angle, alpha);
[U,S,V] = svd(A);
t_orig = original_transmitter_weights(lamda, W, angle, alpha);
error_delta_range = logspace(-15, 0, 15);
%m
[t_tild_list, error_list] = find_t_and_error(A, U, S, V, t_orig, error_delta_range);

figure(1)
scatter(error_delta_range, error_list,'filled','LineWidth',2)
set(gca,'xscale','log','YScale','log')
grid on
hold on
plot(error_delta_range,error_list)
hold off
condition_num_A = cond(A);
disp('kappa value:')
disp(num2str(condition_num_A,5))

%% n - solve with Algo.1

tau_list = [10^-2, 10^-5, 10^-8];
gamma = 30;
error_list_tau = zeros(length(tau_list), length(error_delta_range));
Algo_num = 1; %1 or 2

for t=1:size(tau_list,2)
    tau = tau_list(1,t);
    if Algo_num == 1
        tau_r_or_B_0 = 0.1; %for tau_r
    elseif Algo_num == 2
        tau_r_or_B_0 = 100; %for B_0
    else 
        disp("Wrong Algo_num values")
    end
    [U_LR_hat,S_LR_hat,V_LR_hat, Rank_l] = Improve_accuracy_for_LS_solve(Algo_num, A, gamma, tau, tau_r_or_B_0);

    [t_tild_list ,error_list2] = find_t_and_error(A,U_LR_hat,S_LR_hat,V_LR_hat,t_orig, error_delta_range);
    error_list_tau(t,:) = error_list2;
end
figure(1)
scatter(error_delta_range2,error_list_tau,'filled','LineWidth',2)
set(gca,'xscale','log','YScale','log')
grid on
hold on
%plot(error_delta_range2,error_list_tau)
hold off
condition_num_A = cond(A);
disp('kappa value:')
disp(num2str(condition_num_A,5))
    
    



