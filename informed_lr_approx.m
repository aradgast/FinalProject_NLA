% This function implementing the informed LR approxamtion

function [U, B] = informed_lr_approx(A, gamma)
    tau = 10 ^ -5;
    tau_r = 0.1;
    [Rank_l, n, l] = fast_rank_estimation_LR_approximation(A,tau, tau_r);
    G = randn(size(A, 1), Rank_l + gamma);
    M = A * G;
    [U, S, V] = svd(M);
    U_hat = U(:, 1:Rank_l);
    B = conj(U_hat') * A;
    U = U_hat;
end