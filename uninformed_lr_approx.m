% This function implementing the Un-informed iterative LR approxamtion
function [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma, B_0, tau)
    U_hat = []; B_hat = [];
    l = 1;
    a_norm = norm(A, "fro");
    while true
        
        G = randn(size(A,1), B_0 + gamma) + 1i * randn(size(A,1), B_0 + gamma);
        M = A * G;
        [U, S, V] = svd(M);
        B = conj(U);
        B = B(:, 1:B_0)'*A;
        U_hat = [U_hat U(:, 1:B_0)];
        B_hat = [B_hat; B];
        if norm(B, "fro") <= tau * a_norm
            break;
        else
            disp(norm(B, "fro"));
            A = A - U(:, 1:B_0)*B;
            l = l+1;
        end
    end
    rank_l = l * B_0;
end