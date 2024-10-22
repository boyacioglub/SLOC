% This example calculates 3-sigma outliers for the single axis attitude 
% estimation problem given in Crassidis and Junkins (2012, pp. 165-168) 
% using the recursive formula provided in Tichavsky et al. (1998). The 
% bounds are expected to match Fig. 3.3 (p. 167).  The code is written by
% Burak Boyacioglu in October 2024.

dt = 1;
tf = 60*60;
t = 0:dt:tf;
w = length(t);
N = w(end)-1;

A = [1 -dt; 0 1];
C = [1 0];

% Gyro and Attitude Parameters
sig_u = sqrt(10)*1e-10;
sig_v = sqrt(10)*1e-7;
sig_n = 17*1e-6;
R = sig_n^2;
Q = [sig_v^2*dt+1/3*sig_u^2*dt^3 -1/2*sig_u^2*dt^2;-1/2*sig_u^2*dt^2 sig_u^2*dt];

n = size(A,2); % no. of state variables
p = size(C,1); % no. of outputs

As = zeros(n,n,N+1);
Cs = zeros(p,n,N+1);
Qs = zeros(n,n,N+1);
Rs = zeros(p,p,N+1);

for k=0:N
    As(:,:,k+1) = A;
    Cs(:,:,k+1) = C;
    Qs(:,:,k+1) = Q;
    Rs(:,:,k+1) = R;
end

F_xN = stochConstructabilityGram(As,Cs,Qs,Rs,N);

% Since the system is linear time-invariant, one can calculate the final
% matrix (for a large enough time window) using the following formula. 
F_DARE_final = idare(inv(A),eye(2),C'*inv(R)*C,A'*inv(Q)*A);

P = zeros(n,n,w);
for i=1:w
    P(:,:,i) = inv(F_xN(:,:,i));
end

figure
plot((10:w)/60-1/60,squeeze(3e6*sqrt(P(1,1,10:N+1))),'LineWidth',1.5)
hold on
plot((10:w)/60-1/60,-squeeze(3e6*sqrt(P(1,1,10:N+1))),'LineWidth',1.5)
xlabel('time (minute)')
ylabel('3\sigma bound on x_1 (microradian)')
