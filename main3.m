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
gamma_list = (1:5:150);
M = length(gamma_list);
A_norm = norm(A,'fro');
informed_time = zeros(1,M);
informed_rank = zeros(1,M);
relative_error = zeros(1,M);

tau = 10^-5;
tau_r = 0.1;
mc_num=100;

%get classic truncated SVD algo time
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
%%
figure(1)
scatter(gamma_list,informed_time,'filled','LineWidth',2)
hold on
title('3.k - Time Complexity as a function of \gamma')
grid on
p = polyfit(gamma_list,informed_time,1);
f = polyval(p,gamma_list);
plot(gamma_list,f,"Color",[0,0,1])
plot(gamma_list,simple_tranc_time*ones(1,M),'LineStyle','--',"LineWidth",2)
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
figure(3)
scatter(gamma_list,relative_error,'filled','LineWidth',2)
xlabel("\gamma")
grid on
title('3.k - Relative Error as a function of \gamma')
hold on
plot(gamma_list,relative_error,"Color",[0,0,1])
hold off
%% i
clear
close all
%Parameters
lamda = 1;
W = 4*lamda;
angle = pi/2;
alpha = 1;

%Initialization
test_max = 7;
tao_list = [10^-2,10^-5,10^-8];
T = length(tao_list);
A_truncated_SVD_time = zeros(test_max,T);
A_informed_SVD_time = zeros(test_max,T);
A_mat_memory_store = zeros(test_max,T);
simple_rank = zeros(test_max,T);
informed_rank = zeros(test_max,T);
x_axis = zeros(1,test_max);

tau_r = 0.1;
gamma = 46;
mc_num=100;

for p=1:test_max
    W_test = W*2^p;

    %create A matrix
    A_test = create_steering_mat(lamda, W_test, angle, alpha);
    N_test = size(A_test,1); 
    x_axis(1,p) = N_test;
    
    for t=1:size(tao_list,2)
        tau_test = tao_list(1,t);
        %compute tranked svd time:
        tic;
        [U_tranc,S_tranc,V_tranc] = svd(A_test);
        simple_rank(p,t) = sum(diag(S_tranc)>tau_test);
        A_truncated = U_tranc(:,1:simple_rank)* S_tranc(1:simple_rank,1:simple_rank)* V_tranc(:,1:simple_rank)';
        A_truncated_SVD_time(p,t) = toc;
    
        %compute informed svd algo time:
        tic;
        [U_informed, B_informed, rank_l] = informed_lr_approx(A_test, gamma, tau_test, tau_r);
        A_informed = U_informed*B_informed;
        A_informed_SVD_time(p,t) = toc;
        informed_rank(p,t) = rank_l;       
    end
end
%%
figure()  
hold on
marker = ['x','+','>'];

for t=1:size(tao_list,2)
    scatter(x_axis,A_truncated_SVD_time(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
    scatter(x_axis,A_informed_SVD_time(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
end
set(gca,'xscale','log','yscale','log')
title("3.l - SVD Computation Time as a funtion of N")
xlabel("N")
ylabel('Time [sec]')
legend('\tau = 10^-2, Truncated SVD','\tau = 10^-2, Alg.1','\tau = 10^-5, Truncated SVD', ...
    '\tau = 10^-5, Alg.1','\tau = 10^-8, Truncated SVD''\tau = 10^-8, Alg.1')
grid on
hold off

figure()
hold on
for t=1:size(tao_list,2)
    scatter(x_axis,simple_rank(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
    scatter(x_axis,informed_rank(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
end

title('3.l - Rank as a function of N')
xlabel("N")
ylabel('Rank')
set(gca,'xscale','log','yscale','log')
grid on
legend('\tau = 10^-2, Truncated SVD','\tau = 10^-2, Alg.1','\tau = 10^-5, Truncated SVD', ...
    '\tau = 10^-5, Alg.1','\tau = 10^-8, Truncated SVD','\tau = 10^-8, Alg.1')



