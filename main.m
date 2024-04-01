clc
clear all
close all
%% 1 - Parameters
lamda = 1;
W = 4*lamda;
angle = pi/6;
alpha = 1;

%% 1.a
A = create_steering_mat(lamda, W, angle, alpha);
N = size(A,1); 
figure(1)
imagesc(abs(A))
colorbar;
title("1.a - A Matrix plot")

S = svd(A);
figure(2)
stem(1:N,S,'filled','LineWidth',1)
title("1.a - Singular Values of A")
xlabel("n")
ylabel("Singular Values")
grid on
grid minor

%% 1.b
test_max = 5;
tao_list = [10^-2,10^-5,10^-8];
A_mat_fill_time = zeros(1,test_max);
A_mat_SVD_time = zeros(1,test_max);
A_mat_memory_store = zeros(1,test_max);
rank = zeros(test_max,3);
x_axis = zeros(1,test_max);

for p=1:test_max
    W_test = W^p;
    %create A matrix and compute time
    tic;
    A_test = create_steering_mat(lamda, W_test, angle, alpha);
    A_mat_fill_time(1,p) = toc;

    N_test = size(A_test,1); 
    x_axis(1,p) = N_test;
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
%%
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

%% 1.c
figure(7)
scatter(x_axis,rank,'filled','LineWidth',1)
title('1.c - Rank as a function of N for \tau')
xlabel("N")
ylabel('Rank')
set(gca,'xscale','log','yscale','log')
grid on
legend('\tau = 10^-2','\tau = 10^-5','\tau = 10^-8')

%% 1.d - repeating c for teta=0
angle_d = 0;

test_max = 5;
rank_d = zeros(test_max,3);
x_axis = zeros(1,test_max);

for p=1:test_max
    W_test = W^p;
    %Create A matrix
    A_test = create_steering_mat(lamda, 100*W_test, angle_d, alpha);

    N_test = size(A_test,1); 
    x_axis(1,p) = N_test;
    %Calculate SVD
    S_test = svd(A_test);

    if p==2 || p==4
        figure()
        stem(1:N_test,S_test,'filled','LineWidth',1)
        title("1.d - Singular Values of A for N = " + N_test)
        xlabel("n")
        ylabel("Singular Values")
        grid on
    end

    for t=1:size(tao_list,2)
        rank_d(p,t) = sum(S_test>tao_list(1,t));
    end
end

figure(10)
scatter(x_axis,rank_d,'filled','LineWidth',2)
title('1.d - Rank as a function of N with \theta = 0')
xlabel("N")
ylabel('Rank')
set(gca,'xscale','log')
grid on
legend('\tau = 10^-2','\tau = 10^-5','\tau = 10^-8')

%% 1.e - 
W_e = 64*lamda;
%a - alpha=1, modifying angle
test_max_e = 50;
theta_e = linspace(pi/2,0,test_max_e);
tao_e = 10^-5;

rank_e = zeros(1,test_max_e);

for i=1:test_max_e
    angle_test = theta_e(1,i);
    %Create A matrix
    A_test = create_steering_mat(lamda, W_e, angle_test, alpha);
    N_test = size(A_test,1); 
    %Calculate SVD
    S_test = svd(A_test);

    % if p==2 || p==4
    %     figure()
    %     stem(1:N_test,S_test,'filled','LineWidth',1)
    %     title("1.d - Singular Values of A for N = " + N_test)
    %     xlabel("n")
    %     ylabel("Singular Values")
    % end

    rank_e(1,i) = sum(S_test>tao_e);

end

figure(11)
scatter(theta_e,rank_e,'filled','LineWidth',2)
title('1.e - Rank as a function of \theta ')
xlabel("\theta [rad]")
ylabel('Rank')
xticks(0:pi/6:pi/2);
xticklabels({'0', '\pi/6', '\pi/3', '\pi/2'});
grid on


%b - angle=pi/2, modifying alpha
alpha_e = logspace(-2,2,test_max_e);
rank_e_2 = zeros(1,test_max_e);

for j=1:test_max_e
    alpha_test = alpha_e(1,j);
    %Create A matrix
    A_test = create_steering_mat(lamda, W_e, angle, alpha_test);
    N_test = size(A_test,1); 
    %Calculate SVD
    S_test = svd(A_test);

    % if p==2 || p==4
    %     figure()
    %     stem(1:N_test,S_test,'filled','LineWidth',1)
    %     title("1.d - Singular Values of A for N = " + N_test)
    %     xlabel("n")
    %     ylabel("Singular Values")
    % end

    rank_e_2(1,j) = sum(S_test>tao_e);
end

figure(12)
scatter(alpha_e,rank_e_2,'filled','LineWidth',2)
title('1.e - Rank as a function of \alpha ')
xlabel("\alpha")
ylabel('Rank')
set(gca,'xscale','log')
grid on
