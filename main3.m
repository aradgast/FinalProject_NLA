clc
clear
close all

%% 3 - Parameters
lamda = 1;
W = 64*lamda;
angle = pi/2;
alpha = 1;
A = create_steering_mat(lamda, W, angle, alpha);
%%
B_0 = 5;
tic;
S = svd(A);
svd_time = toc;
mc_num=100;

svd_rank = zeros(1, 3);
informed_time = zeros(2, 3);
informed_rank = zeros(2, 3);

idx_r = 0;
idx = 0;
gamma = 100;

for tau_r = [0.1, 0.01]
    idx_r = idx_r + 1;
    for tau=[10^-2, 10^-5, 10^-8]
        idx = idx + 1;
        informed_rank_mc = zeros(1,mc_num);
        for mc=1:mc_num
            tic;
            [U_hat, B_hat, rank_l] = informed_lr_approx(A, gamma, tau, tau_r);
            informed_time(idx_r, idx) = informed_time(idx_r, idx)+toc/mc_num;
            informed_rank_mc(1,mc) = rank_l;
        end
        informed_rank(idx_r, idx) = median(informed_rank_mc);
        if idx_r == 1
            svd_rank(1, idx) = sum(S > tau);
        end
    end
    idx = 0;
end

%% Gamma dependent time complexity and accuracy
gamma_list = (1:20:600);
M = length(gamma_list);
A_norm = norm(A,'fro');
informed_time = zeros(1,M);
simple_tranc_time = zeros(1,M);

informed_rank = zeros(1,M);
relative_error = zeros(1,M);

tau = 10^-5;
tau_r = 0.1;
mc_num=100;

tic;
[U,S,V] = svd(A);
simple_rank = sum(diag(S)>tau);
A_truncated = U(:,1:simple_rank)* S(1:simple_rank,1:simple_rank)* V(:,1:simple_rank)';
simple_tranc_time = toc;

for g=1:length(gamma_list)
    gamma_test = gamma_list(1,g);
    relative_error_mc = zeros(1,mc_num);

    for mc=1:mc_num
        %get algo time
        tic;
        [U_hat, B_hat, rank_l] = informed_lr_approx(A, gamma_test, tau, tau_r);
        A_informed = U_hat*B_hat;
        end_informed_algo_time = toc;
        informed_time(1, g) = informed_time(1, g)+end_informed_algo_time/mc_num;
        
        %relative error
        relative_error(1,g) = relative_error(1,g) + (norm(A_truncated-A_informed,'fro')/norm(A_truncated,'fro'))/mc_num;
    end

end

figure(3)
scatter(gamma_list,relative_error,'filled','LineWidth',2)
xlabel("\gamma")
grid on
title('3.k - Relative Error as a function of \gamma')
hold on
p = polyfit(gamma_list,relative_error,1);
f = polyval(p,gamma_list);
plot(gamma_list,f,"Color",[0,0,1])
plot(gamma_list,relative_error,"Color",[0,0,1])
legend('Error Values', 'Trend Line')
hold off

figure(1)
scatter(gamma_list,informed_time,'filled','LineWidth',2)
hold on
title('3.k - Time Complexity as a function of \gamma')
grid on
p = polyfit(gamma_list,informed_time,1);
f = polyval(p,gamma_list);
plot(gamma_list,f,"Color",[0,0,1])
plot(gamma_list,simple_tranc_time*ones(1,M),'LineStyle','--')
hold off
xlabel("\gamma")
ylabel("Time [sec]")
legend('Fast LR Approx Time Points', 'Trend Line','Straightforward SVD')
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

%% i
test_max = 10;
tao_list = [10^-2,10^-5,10^-8];
A_mat_fill_time = zeros(1,test_max);
A_mat_SVD_time = zeros(1,test_max);
A_mat_memory_store = zeros(1,test_max);
rank = zeros(test_max,3);
x_axis = zeros(1,test_max);

tau_r = 0.1;
gamma = 100;
mc_num=100;

for p=1:test_max
    W_test = W*2^p;
    %create A matrix
    A_test = create_steering_mat(lamda, W_test, angle, alpha);

    N_test = size(A_test,1); 
    x_axis(1,p) = N_test;
    
    %compute clasic SVD:

    %Compute memory storage in bytes
    A_mat_memory_store(1,p) = 8*2*N_test^2;

    tic;
    S_test = svd(A_test);
    A_mat_SVD_time(1,p) = toc;
    if p==2 || p==4
        figure()
        stem(1:N_test,S_test,'filled','LineWidth',1)
        title("1.c - Singular Values of A for N = " + N_test)
        xlabel("n")
        ylabel("Singular Values")
        grid on
    end

    for t=1:size(tao_list,2)
        rank(p,t) = sum(S_test>tao_list(1,t));
    end
end

figure(5)
scatter(x_axis,A_mat_fill_time,'filled','LineWidth',2)
hold on
scatter(x_axis,A_mat_SVD_time,'filled','LineWidth',2)
set(gca,'xscale','log','yscale','log')
title("1.b - Computation Time as a funtion of N")
xlabel("N")
ylabel('Time [sec]')
legend('A fill', 'SVD')
grid on
hold off

figure(6)
scatter(x_axis,A_mat_memory_store,'filled','LineWidth',2)
set(gca,'xscale','log','yscale','log')
title("1.b - Memory Storage as a funtion of N")
xlabel("N")
ylabel('Memory Storage [bytes]')
grid on

figure(7)
scatter(x_axis,rank,'filled','LineWidth',1)
title('1.c - Rank as a function of N for \tau')
xlabel("N")
ylabel('Rank')
set(gca,'xscale','log','yscale','log')
grid on
legend('\tau = 10^-2','\tau = 10^-5','\tau = 10^-8')


