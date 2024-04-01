clc
clear all
close all
%% 1 - Parameters
lamda = 1;
W = 64*lamda;
angle = pi/2;
alpha = 1;
tao = 10^-5;
%% 2.h
A = create_steering_mat(lamda, W, angle, alpha);
tau_r = [10^-1, 5*10^-2, 10^-2];
tau = [10^-2, 10^-5, 10^-8];
A_Ranks = zeros(1,length(tau));
N_values = size(A,1)*ones(1,length(tau));
Ranks = zeros(length(tau_r),length(tau));
n_values = zeros(length(tau_r),length(tau));
iter = zeros(length(tau_r),length(tau));

for b=1:length(tau)
    for a=1:length(tau_r)
        [Ranks(a,b), n_values(a,b),iter(a,b)] = fast_rank_estimation_LR_approximation(A,tau(b), tau_r(a));
    end
    A_Ranks(1,b) = sum(svd(A)>tau(b));
end 

% Create a table
T1 = table(tau', A_Ranks', N_values', Ranks(1,:)', n_values(1,:)',iter(1,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
T2 = table(tau', A_Ranks', N_values', Ranks(2,:)', n_values(2,:)',iter(2,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
T3 = table(tau', A_Ranks', N_values', Ranks(3,:)', n_values(3,:)',iter(3,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
%T2.NewColumn = tau_r(1,2);
% Display the table
disp(T1);
disp(T2);
disp(T3);

%% 2.i.1
W_2 = 128*lamda;
A = create_steering_mat(lamda, W_2, angle, alpha);
tau_r = [10^-1, 5*10^-2, 10^-2];
tau = [10^-2, 10^-5, 10^-8];
A_Ranks = zeros(1,length(tau));
N_values = size(A,1)*ones(1,length(tau));
Ranks = zeros(length(tau_r),length(tau));
n_values = zeros(length(tau_r),length(tau));
iter = zeros(length(tau_r),length(tau));

for b=1:length(tau)
    for a=1:length(tau_r)
        [Ranks(a,b), n_values(a,b),iter(a,b)] = fast_rank_estimation_LR_approximation(A,tau(b), tau_r(a));
    end
    A_Ranks(1,b) = sum(svd(A)>tau(b));
end 

% Create a table
T1 = table(tau', A_Ranks', N_values', Ranks(1,:)', n_values(1,:)',iter(1,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
T2 = table(tau', A_Ranks', N_values', Ranks(2,:)', n_values(2,:)',iter(2,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
T3 = table(tau', A_Ranks', N_values', Ranks(3,:)', n_values(3,:)',iter(3,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
%T2.NewColumn = tau_r(1,2);
% Display the table
disp(T1);
disp(T2);
disp(T3);

%% 2.i.2
angle_2 = pi/2;
A = create_steering_mat(lamda, W, angle_2, alpha);
tau_r = [10^-1, 5*10^-2, 10^-2];
tau = [10^-2, 10^-5, 10^-8];
A_Ranks = zeros(1,length(tau));
N_values = size(A,1)*ones(1,length(tau));
Ranks = zeros(length(tau_r),length(tau));
n_values = zeros(length(tau_r),length(tau));
iter = zeros(length(tau_r),length(tau));

for b=1:length(tau)
    for a=1:length(tau_r)
        [Ranks(a,b), n_values(a,b),iter(a,b)] = fast_rank_estimation_LR_approximation(A,tau(b), tau_r(a));
    end
    A_Ranks(1,b) = sum(svd(A)>tau(b));
end 

% Create a table
T1 = table(tau', A_Ranks', N_values', Ranks(1,:)', n_values(1,:)',iter(1,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
T2 = table(tau', A_Ranks', N_values', Ranks(2,:)', n_values(2,:)',iter(2,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
T3 = table(tau', A_Ranks', N_values', Ranks(3,:)', n_values(3,:)',iter(3,:)', ...
    'VariableNames', {'Tau ', 'A Rank', 'N', 'RankFound', 'n', 'Iterantions'});
%T2.NewColumn = tau_r(1,2);
% Display the table
disp(T1);
disp(T2);
disp(T3);