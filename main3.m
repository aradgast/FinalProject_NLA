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


%%
uninformed_time = zeros(1, 3);
uninformed_rank = zeros(1, 3);