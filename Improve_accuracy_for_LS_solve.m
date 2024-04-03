function [U_LR_hat,S_LR_hat,V_LR_hat, Rank_l] = Improve_accuracy_for_LS_solve(num, A, gamma, tau, tau_r_or_B_0)
%proposed to improve the accuracy by solving the system, in the least squares sense, using a
%truncated SVD of A, with arbitrary tau
%num is to decide which SVD estimation to use. It can be equal to only 1 or 2

if num==1 %use informed algorithm
    tau_r = tau_r_or_B_0;
    [U_hat, B_hat, Rank_l] = informed_lr_approx(A, gamma, tau, tau_r); A_approx = U_hat*B_hat;

elseif num==2 %use uninformed algorithm
    B_0 = tau_r_or_B_0;
    [U_hat, B_hat, Rank_l] = uninformed_lr_approx(A, gamma, B_0, tau);

else %use informed algorithm
    disp("Wrong input- please choose a 1 for Algo.1 or 2 for Algo.2")
    return
end

[U_B, S_B, V_B] = svd(B_hat);
U_LR_hat = U_hat*U_B;
S_LR_hat = S_B(1:Rank_l,1:Rank_l);
V_LR_hat = V_B;

end