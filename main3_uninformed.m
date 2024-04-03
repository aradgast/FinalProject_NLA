clc
clear all
close all

%% 3 - Parameters
lamda = 1;
W = 64*lamda;
angle = pi/2;
alpha = 1;
A = create_steering_mat(lamda, W, angle, alpha);
%%
B_0 = 15;
tic;
S = svd(A);
svd_time = toc;
svd_rank = zeros(1, 3);
un_informed_time = zeros(2, 3);
uninformed_rank = zeros(2, 3);

idx = 0;
idx_g = 0;
monte_carlo_num = 100;
for gamma = [10, 50]
    idx_g = idx_g + 1;
    for tau=[10^-2, 10^-5, 10^-8]
        idx = idx + 1;
        un_informed_rank_tmp = zeros(1, monte_carlo_num);
        for mc=1:monte_carlo_num
            tic;
            [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma, B_0, tau);
            un_informed_time(idx_g, idx) = un_informed_time(idx_g, idx) + toc / monte_carlo_num;
            un_informed_rank_tmp(1, mc) = rank_l;
        end
        uninformed_rank(idx_g, idx) = median(un_informed_rank_tmp);
        if idx_g == 1
            svd_rank(1, idx) = sum(S > tau);
        end
    end
    idx = 0;
end


%% Gamma dependent time complexity and accuracy
b_list = (1:1:80);
M = length(b_list);
A_norm = norm(A,'fro');
uninformed_time = zeros(1,M);
simple_tranc_time = zeros(1,M);

uninformed_rank = zeros(1,M);
relative_error = zeros(1,M);

tau = 10^-5;
gamma = 40;
mc_num=2;

tic;
[U,S,V] = svd(A);
simple_rank = sum(diag(S) > tau);
A_truncated = U(:,1:simple_rank)* S(1:simple_rank,1:simple_rank)* V(:,1:simple_rank)';
simple_tranc_time = toc;
A_truncated_norm = norm(A_truncated,'fro');

for idx_b=1:length(b_list)
    b_test = b_list(1,idx_b);

    for mc=1:mc_num
        %get algo time
        tic;
        [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma, b_test, tau);
        A_uninformed = U_hat * B_hat;
        end_uninformed_algo_time = toc;
        uninformed_time(1, idx_b) = uninformed_time(1, idx_b) + end_uninformed_algo_time / mc_num;
        
        %relative error
        relative_error(1, idx_b) = relative_error(1, idx_b) + (norm(A_truncated-A_uninformed,'fro') / A_truncated_norm) / mc_num;
    end

end

figure(3)
scatter(b_list,relative_error,'filled','LineWidth',2)
xlabel("B_{0}")
grid on
title('3.k - Relative Error as a function of B_{0}')
hold on
p = polyfit(b_list,relative_error,1);
f = polyval(p,b_list);
plot(b_list,f,"Color",[0,0,1])
plot(b_list,relative_error,"Color",[0,0,1])
legend('Error Values', 'Trend Line')
hold off

figure(1)
scatter(b_list,uninformed_time,'filled','LineWidth',2)
hold on
title('3.k - Time Complexity as a function of B_{0}')
grid on
p = polyfit(b_list,uninformed_time,1);
f = polyval(p,b_list);
plot(b_list,f,"Color",[0,0,1])
plot(b_list,simple_tranc_time*ones(1,M),'LineStyle','--')
hold off
xlabel("B_{0}")
ylabel("Time [sec]")
legend('Fast LR Approx Time Points', 'Trend Line','Straightforward SVD')
% figure(2)
% scatter(gamma_list,fast_rank_estimation,'filled','LineWidth',2)
% hold on
% scatter(gamma_list,informed_rank,'filled','LineWidth',2)
% scatter(gamma_list,Orig_rank.*ones(1,M),'LineWidth',3)
% title('Rank estimation as a function of \gamma')
% xlabel("\gamma")
% grid on
% legend('Fast Rank Estimation algorithm', 'Fast A estimation', 'Real rank')
% hold off





