function t = original_transmitter_weights(lamda, W, angle, alpha)
k = 2*pi/lamda;
k0 = k/4;
delta = lamda/10;
N = W/delta;

r_t1 = zeros(2,N); %transmitter 1 (x,y) compoinent
r_t = zeros(2,N); %transmitter (x,y) compoinents 

r_t(1,:) = (-1)*linspace(alpha*W/2,alpha*W/2+N*delta*cos(angle),N);
r_t(2,:) = linspace(0,N*delta*sin(angle),N);
r_t1 = r_t(:,1);

diff_x = r_t(1,:)-r_t1(1,1); %1xN
diff_y = r_t(2,:)-r_t1(2,1); %1xN

diff = sqrt(diff_x.^2 + diff_y.^2); %should be 1xN

t = (exp(1i*k0*diff))'; %should be Nx1
