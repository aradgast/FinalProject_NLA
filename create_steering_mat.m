function A = create_steering_mat(lamda, W, angle, alpha)
delta = lamda/10;
N = W/delta;
A = zeros(N,N);
%creat reciever transmitter
r_r = zeros(2,N); %reciever and transmitter are symetric. each with (x,y) compoinents 
r_r(1,:) = linspace(alpha*W/2,alpha*W/2+N*delta*cos(angle),N);
r_r(2,:) = linspace(0,N*delta*sin(angle),N);

r_t = r_r;
r_t(1,:) = r_t(1,:)*(-1);

diff_x = r_r(1,:)-r_t(1,:)';
diff_y = r_r(2,:)-r_t(2,:)';
diff = sqrt(diff_x.^2 + diff_y.^2); %should be NxN matrix

k=2*pi/lamda;
A = exp(-1i*k.*diff)./sqrt(diff); %should be NxN matrix

end