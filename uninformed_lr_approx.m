% This function implementing the Un-informed iterative LR approxamtion
function [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma, B_0)
    U_hat = []; B_hat = [];
    a_norm = norm(A, "fro");
    l = 1;
    while true
        G = randn(size(A,1), B_0 + gamma);
        M = A * G;
        [U, S, V] = svd(M);
        B = conj(U)(:, 1:B_0);
        U_hat = [U_hat U(:, 1:B_0)];
        B_hat = [B_hat B];
        if norm(B, "fro") <= tau * a_norm
            break;
        else
            A = A - U(:, 1:B_0)*B;
            l = l+1;
        end
    end
    rank_l = l * B_0;
end