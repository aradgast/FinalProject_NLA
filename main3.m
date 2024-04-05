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
tic;
S = svd(A);
svd_time = toc;
mc_num=99;

tau_r_l = [0.1, 0.01];
tau_l = [10^-2, 10^-5, 10^-8];

svd_rank = zeros(1, length(tau_l));
informed_time = zeros(length(tau_r_l), length(tau_l));
informed_rank = zeros(length(tau_r_l), length(tau_l));

idx_r = 0;
idx = 0;
gamma = 10;

for tau_r = tau_r_l
    idx_r = idx_r + 1;
    for tau= tau_l
        idx = idx + 1;
        informed_rank_mc = zeros(1,mc_num);
        for mc=1:mc_num
            tic;
            [U_hat, B_hat, rank_l] = informed_lr_approx(A, gamma, tau, tau_r);
            end_time = toc;
            informed_time(idx_r, idx) = informed_time(idx_r, idx)+end_time/mc_num;
            informed_rank_mc(1,mc) = rank_l;
        end
        informed_rank(idx_r, idx) = median(informed_rank_mc);
        if idx_r == 1
            svd_rank(1, idx) = sum(S > tau);
        end
    end
    idx = 0;
end

%% k - Gamma dependent time complexity and accuracy
gamma_list = (1:5:150);
M = length(gamma_list);
A_norm = norm(A,'fro');
informed_time = zeros(1,M);
relative_error = zeros(1,M);

tau = 10^-5;
tau_r = 0.1;
mc_num=99;

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
legend('Fast LR Approx Time Points', 'Trend Line','Straightforward SVD', 'Location','west')
f = gcf;
exportgraphics(f,'3.k - Time Complexity as a function of gamma.jpg','Resolution',150)

figure(2)
scatter(gamma_list,relative_error,'filled','LineWidth',2)
xlabel("\gamma")
grid on
title('3.k - Relative Error as a function of \gamma')
hold on
plot(gamma_list,relative_error,"Color",[0,0,1])
hold off
f = gcf;
exportgraphics(f,'3.k - Relative Error as a function of gamma.jpg','Resolution',150)
%% i
clear

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
mc_num=49;

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
close all
figure(3)  
hold on
marker = ['x','+','>','x','+','>'];
SVDfor_plot = [A_truncated_SVD_time, A_informed_SVD_time];
for t=1:size(SVDfor_plot,2)
    if t <4
    scatter(x_axis,SVDfor_plot(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
    else
    scatter(x_axis,SVDfor_plot(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
    end
end
set(gca,'xscale','log','yscale','log')
title("3.l - Computation Time as a funtion of N")
xlabel("N")
ylabel('Time [sec]')
legend('\tau = 10^-2, Truncated SVD','\tau = 10^-5, Truncated SVD','\tau = 10^-8, Truncated SVD' ...
    ,'\tau = 10^-2, Alg.1','\tau = 10^-5, Alg.1','\tau = 10^-8, Alg.1')
grid on
hold off
f = gcf;
exportgraphics(f,'3.l - Computation Time as a funtion of N.jpg','Resolution',150)

%Use this to zoom in to specific N values:
% zoom_N = 3;
% figure(4)  
% hold on
% marker = ['x','+','>','x','+','>'];
% SVDfor_plot = [A_truncated_SVD_time, A_informed_SVD_time];
% for t=1:size(SVDfor_plot,2)
%     if t <4
%     scatter(x_axis(1,zoom_N:zoom_N),SVDfor_plot(zoom_N:zoom_N,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
%     else
%     scatter(x_axis(1,zoom_N:zoom_N),SVDfor_plot(zoom_N:zoom_N,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
%     end
% end
% hold off
% set(gca,'xscale','log','yscale','log')
% title("3.l - SVD Computation Time as a funtion of N")
% xlabel("N")
% ylabel('Time [sec]')
% legend('\tau = 10^-2, Truncated SVD','\tau = 10^-5, Truncated SVD','\tau = 10^-8, Truncated SVD' ...
%     ,'\tau = 10^-2, Alg.1','\tau = 10^-5, Alg.1','\tau = 10^-8, Alg.1')
% grid on


figure(5)
hold on
rank_for_plot = [simple_rank,informed_rank];
for t=1:size(rank_for_plot,2)
    if t <4
    scatter(x_axis,rank_for_plot(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0])
    else
    scatter(x_axis,rank_for_plot(:,t)',marker(t),'LineWidth',1.5,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1])
    end
end

title('3.l - Rank as a function of N')
xlabel("N")
ylabel('Rank')
set(gca,'xscale','log','yscale','log')
grid on
legend('\tau = 10^-2, Truncated SVD','\tau = 10^-5, Truncated SVD','\tau = 10^-8, Truncated SVD' ...
    ,'\tau = 10^-2, Alg.1','\tau = 10^-5, Alg.1','\tau = 10^-8, Alg.1')
f = gcf;
exportgraphics(f,'3.l - Rank as a function of N.jpg','Resolution',150)


