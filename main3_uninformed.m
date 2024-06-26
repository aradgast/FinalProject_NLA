clc
clear
close all

% 3 - Parameters
lamda = 1;
W = 64*lamda;
angle = pi/2;
alpha = 1;
A = create_steering_mat(lamda, W, angle, alpha);
%% Understanding the algorithm with different parameters
B_0 = 15;
tic;
S = svd(A);
svd_time = toc;
gamma_list = [10, 150];
tau_list = [10^-1, 10^-2, 10^-5];
svd_rank = zeros(1, 3);
un_informed_time = zeros(length(gamma_list), length(tau_list));
uninformed_rank = zeros(length(gamma_list), length(tau_list));

idx = 0;
idx_g = 0;
monte_carlo_num = 99;
for gamma = gamma_list
    idx_g = idx_g + 1;
    for tau=tau_list
        idx = idx + 1;
        un_informed_rank_tmp = zeros(1, monte_carlo_num);
        for mc=1:monte_carlo_num
            tic;
            [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma, B_0, tau);
            un_informed_time(idx_g, idx) = un_informed_time(idx_g, idx) + toc / monte_carlo_num;
            un_informed_rank_tmp(1, mc) = rank_l;
        end
        uninformed_rank(idx_g, idx) = median(un_informed_rank_tmp);
        if idx_g == 1
            svd_rank(1, idx) = sum(S > tau);
        end
    end
    idx = 0;
end

for idx_g = 1:length(gamma_list)
    for idx_tau = 1:length(tau_list)
        gamma_val = gamma_list(idx_g);
        tau_val = tau_list(idx_tau);
        
        % Extract results for current combination of gamma and tau
        un_informed_time_val = un_informed_time(idx_g, idx_tau);
        uninformed_rank_val = uninformed_rank(idx_g, idx_tau);
        svd_rank_val = svd_rank(idx_tau);
        
        % Create table for current combination
        results_table = table(gamma_val, tau_val, un_informed_time_val, uninformed_rank_val, svd_time, svd_rank_val, ...
            'VariableNames', {'Gamma', 'Tau', 'Uninformed_Time', 'Uninformed_Rank', 'SVD_Time', 'SVD_Rank'});
        
        % Display the table
        disp(results_table);
    end
end

%% k - Gamma dependent time complexity and accuracy
b_list = (1:1:15);
M = length(b_list);
tau = 10^-1;
gamma_list = [1,10,100];
G = length(gamma_list);
mc_num=99;

A_norm = norm(A,'fro');
uninformed_time = zeros(G,M);
relative_error = zeros(G,M);

%get classic truncated SVD algo time
tic;
[U,S,V] = svd(A);
simple_rank = sum(diag(S) > tau);
A_truncated = U(:,1:simple_rank)* S(1:simple_rank,1:simple_rank)* V(:,1:simple_rank)';
simple_tranc_time = toc;

A_truncated_norm = norm(A_truncated,'fro');

%get uninformed interative LR approx time and rank
for g=1:G
    gamma_test = gamma_list(1,g);
    for idx_b=1:length(b_list)
        b_test = b_list(1,idx_b);
    
        for mc=1:mc_num
            %get algo time
            tic;
            [U_hat, B_hat, rank_l] = uninformed_lr_approx(A, gamma_test, b_test, tau);
            A_uninformed = U_hat * B_hat;
            end_uninformed_algo_time = toc;
            uninformed_time(g, idx_b) = uninformed_time(g, idx_b) + end_uninformed_algo_time / mc_num;
            
            %relative error
            relative_error(g, idx_b) = relative_error(g, idx_b) + (norm(A_truncated-A_uninformed,'fro') / A_truncated_norm) / mc_num;
        end
    
    end
end
%%
figure(1)
plot(b_list,uninformed_time,'-*','LineWidth',2)
hold on
plot(b_list,simple_tranc_time*ones(1,M),'LineStyle','--')
hold off
title('3.k - Fast LR Approx Time Complexity as a function of B_{0}')
grid on
% p = polyfit(b_list,uninformed_time,1);
% f = polyval(p,b_list);
% plot(b_list,f,"Color",[0,0,1])
xlabel("B_{0}")
ylabel("Time [sec]")
% legend 
legend_str = cell(1, numel(gamma_list));
for i = 1:numel(gamma_list)
    legend_str{i} = sprintf('\\gamma = %d', gamma_list(i));
end
legend(legend_str);
set(gca,'yscale','log')
f = gcf;
exportgraphics(f,'3.k - Uninformed Fast LR Approx Time Complexity as a function of B_0.jpg','Resolution',150)

figure(2)
plot(b_list,relative_error,'-*','LineWidth',2)
xlabel("B_{0}")
ylabel("E_r")
grid on
title('3.k - Relative Error as a function of B_{0}')
% legend 
legend_str = cell(1, numel(gamma_list));
for i = 1:numel(gamma_list)
    legend_str{i} = sprintf('\\gamma = %d', gamma_list(i));
end
legend(legend_str);
set(gca,'yscale','log')
f = gcf;
exportgraphics(f,'3.k - Uninformed Relative Error as a function of B_0.jpg','Resolution',150)

%% i
clear

%Parameters
lamda = 1;
W = 4*lamda;
angle = pi/2;
alpha = 1;

%Initialization
test_max = 7;
tau_list = [10^-2,10^-5,10^-8];
T = length(tau_list);
A_truncated_SVD_time = zeros(test_max,T);
A_uninformed_SVD_time = zeros(test_max,T);
simple_rank = zeros(test_max,T);
uninformed_rank = zeros(test_max,T);
x_axis = zeros(1,test_max);

tau_r = 0.1;
gamma = 46;
B_0 = 10;
mc_num=100;

for p=1:test_max
    W_test = W*2^p;

    %create A matrix
    A_test = create_steering_mat(lamda, W_test, angle, alpha);
    N_test = size(A_test,1); 
    x_axis(1,p) = N_test;
    
    for t=1:size(tau_list,2)
        tau_test = tau_list(1,t);
        %compute tranked svd time:
        tic;
        [U_tranc,S_tranc,V_tranc] = svd(A_test);
        simple_rank(p,t) = sum(diag(S_tranc)>tau_test);
        A_truncated = U_tranc(:,1:simple_rank)* S_tranc(1:simple_rank,1:simple_rank)* V_tranc(:,1:simple_rank)';
        A_truncated_SVD_time(p,t) = toc;
    
        %compute informed svd algo time:
        tic;
        [U_uninformed, B_uninformed, rank_l] = uninformed_lr_approx(A_test, gamma, B_0,tau_test);
        A_uninformed = U_uninformed*B_uninformed;
        A_uninformed_SVD_time(p,t) = toc;
        uninformed_rank(p,t) = rank_l;       
    end
end

figure(3)  
hold on
marker = ['x','+','>'];

for t=1:size(tau_list,2)
    scatter(x_axis,A_truncated_SVD_time(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
    scatter(x_axis,A_uninformed_SVD_time(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
end
set(gca,'xscale','log','yscale','log')
title("3.l - SVD Computation Time as a funtion of N")
xlabel("N")
ylabel('Time [sec]')
legend('\tau = 10^-2, Truncated SVD','\tau = 10^-2, Alg.1','\tau = 10^-5, Truncated SVD', ...
    '\tau = 10^-5, Alg.1','\tau = 10^-8, Truncated SVD''\tau = 10^-8, Alg.1')
grid on
hold off

figure(4)
hold on
for t=1:size(tau_list,2)
    scatter(x_axis,simple_rank(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
    scatter(x_axis,uninformed_rank(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
end

title('3.l - Rank as a function of N')
xlabel("N")
ylabel('Rank')
set(gca,'xscale','log','yscale','log')
grid on
legend('\tau = 10^-2, Truncated SVD','\tau = 10^-2, Alg.1','\tau = 10^-5, Truncated SVD', ...
    '\tau = 10^-5, Alg.1','\tau = 10^-8, Truncated SVD','\tau = 10^-8, Alg.1')








