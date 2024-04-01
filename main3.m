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
gamma = 5;
B_0 = 5;
tic;
[U, S, V] = svd(A);
svd_time = toc;
svd_rank = zeros(1, 3);
informed_time = zeros(2, 3);
informed_rank = zeros(2, 3);

idx_r = 0;
idx = 0;
gamma = 10;
for tau_r = [0.1, 0.01]
    idx_r = idx_r + 1;
    for tau=[10^-2, 10^-5, 10^-8]
        idx = idx + 1;
        tic;
        [U_hat, B_hat, rank_l] = informed_lr_approx(A, gamma, tau, tau_r);
        informed_time(idx_r, idx) = toc;
        informed_rank(idx_r, idx) = rank_l;
        if idx_r == 1
            svd_rank(1, idx) = sum(diag(S) > tau);
        end
    end
    idx = 0;
end


%% Gamma dependent time complexity and accuracy
gamma_list = (1:10:600);
M = length(gamma_list);
N = size(A,1);
A_norm = norm(A,'fro');
informed_time = zeros(1,M);
fast_rank_estimation = zeros(1,M);
informed_rank = zeros(1,M);
relative_error = zeros(1,M);

tau = 10^-5;
tau_r = 0.1;

for k=1:length(gamma_list)
    gamma_test = gamma_list(1,k);
    tic;
    [U_hat, B_hat, rank_l] = informed_lr_approx(A, gamma_test, tau, tau_r);
    informed_time(1, k) = toc;
    fast_rank_estimation(1, k) = rank_l;

    A_hat = U_hat*B_hat;
    relative_error(1,k) = norm(A-A_hat,'fro')/A_norm;

    rank_hat = sum(svd(A_hat)>tau);
    informed_rank(1, k) = rank_hat;   
end
[U, S, V] = svd(A);
Orig_rank = sum(diag(S)>tau);
figure(1)
scatter(gamma_list,informed_time,'filled','LineWidth',2)
hold on
title('Time Complexity as a function og \gamma')
grid on
p = polyfit(gamma_list,informed_time,1);
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




