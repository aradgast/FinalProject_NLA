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
un_informed_rank = zeros(2, 3);

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
        un_informed_rank(idx_g, idx) = median(un_informed_rank_tmp);
        if idx_g == 1
            svd_rank(1, idx) = sum(S > tau);
        end
    end
    idx = 0;
end


%% Gamma dependent time complexity and accuracy
gamma_list = (1:20:600);
B_0 = 10;
tau = 10 ^ -5;
for idx_g=1:length(gamma_list)
    un_informed_rank_tmp = zeros(1, monte_carlo_num);
    for mc=1:monte_carlo_num
        tic;
        [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma, B_0, tau);
        un_informed_time(idx_g, idx) = un_informed_time(idx_g, idx) + toc / monte_carlo_num;
        un_informed_rank_tmp(1, mc) = rank_l;
    end
    un_informed_rank(idx_g, idx) = median(un_informed_rank_tmp);
    if idx_g == 1
        svd_rank(1, idx) = sum(S > tau);
    end
end


figure(1)
scatter(gamma_list,un_informed_time,'filled','LineWidth',2)
hold on
title('Time Complexity as a function og \gamma')
grid on
p = polyfit(gamma_list,un_informed_time,1);
f = polyval(p,gamma_list);
plot(gamma_list,f)
hold off
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
figure(3)
scatter(gamma_list,relative_error,'filled','LineWidth',2)
xlabel("\gamma")
grid on
title('Relative error for A as a function of \gamma')




