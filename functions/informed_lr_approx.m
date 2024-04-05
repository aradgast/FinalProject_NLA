% This function implementing the informed LR approxamtion

function [U_hat, B_hat, Rank_l] = informed_lr_approx(A, gamma, tau, tau_r)
    [Rank_l, ~, ~] = fast_rank_estimation_LR_approximation(A,tau, tau_r);
    G = (randn(size(A, 1), Rank_l + gamma) + 1i * randn(size(A, 1), Rank_l + gamma))/sqrt(2);
    M = A * G;
    [U, ~, ~] = svd(M);
    U_hat = U(:, 1:Rank_l);
    B_hat = U_hat' * A;
end