function [Rank_l,n,l] = fast_rank_estimation_LR_approximation(A,tau, tau_r)
%fast estimation of the rank ℛ(tau) of matrix A for a LR approximation 
%relative error threshold tau.
l=1;
n=10;
N = size(A,1);
Rank_l = sum(svd(A)>tau); %should not effect anything

while n <= N
    Rank_l_prev = Rank_l; 
    %create index vectore
    i_row = randperm(N,n);
    i_col = randperm(N,n);
    
    A_submat_l = A(i_row,i_col);
    S = svd(A_submat_l);
    Rank_l = sum(S>tau);
    % figure()
    % stem(1:length(S),S,'filled','LineWidth',1)
    % xlabel("n")
    % ylabel("Singular Values")
    % grid on
    if l ~= 1
        error_Rank = (Rank_l-Rank_l_prev)/Rank_l;
        if error_Rank< tau_r
            fprintf('error taget reached with Tau = %u and Tau(r) = %u \n', tau, tau_r);
            break
        end 
    end 
    l=l+1;
    n=2*n;
 end 

 if n > N
     n = n/2;
     fprintf('error taget NOT reached with Tau = %u and Tau(r) = %u \n', tau, tau_r);
 end
