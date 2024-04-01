% This function implementing the informed LR approxamtion

function [U, B, Rank_l] = informed_lr_approx(A, gamma, tau, tau_r)
    [Rank_l, n, l] = fast_rank_estimation_LR_approximation(A,tau, tau_r);
    G = randn(size(A, 1), Rank_l + gamma) + 1i * randn(size(A, 1), Rank_l + gamma);
    M = A * G;
    [U, S, V] = svd(M);
    U_hat = U(:, 1:Rank_l);
    B = conj(U_hat') * A;
    U = U_hat;
end