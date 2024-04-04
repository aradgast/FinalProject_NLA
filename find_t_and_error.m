function [t_tild_list, error_list] = find_t_and_error(A,U,S,V,t_orig,error_delta_range)

r = A * t_orig; %Nx1
N = size(A, 1);
% m
error_list = zeros(1, length(error_delta_range));
t_tild_list = zeros(size(r,1), length(error_delta_range));

for d = 1:length(error_delta_range)
    error_delta = error_delta_range(1, d);
    teta = unifrnd(0, 2 * pi, N, 1);
    r_plagued = r + error_delta * (abs(r).*exp(1i * teta));
    S_inv = S^(-1);
    t_plagued = V*S_inv*U'*r_plagued; %Nx1
    t_tild_list(:,d) = t_plagued;
    error_list(1,d) = norm(t_orig-t_plagued,2);
end