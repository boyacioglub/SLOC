function [F_xw] = stochConstructabilityGram(As,Cs,Qs,Rs,N)
% This function generates the w-step constructability Gramians from w=1 to
% w=N+1 for a discrete-time linear time-varying systems with Gaussian noise
% using the recursive formula given in Tichavsky et al. (1998).

% INPUTS
%   As is the collection n-by-n system matrices
%   Cs is the collection p-by-n output matrices
%   Qs is the collection n-by-n process noise covariance matrices
%   Rs is the collection n-by-n measurement noise covariance matrices
%   w is the number of time steps
% OUTPUTS
%   F_x0 is the stochastic constructability Gramian
    n = size(As,1); % number of states
    F_xw = zeros(n,n,N+1);
    F_xw(:,:,1) = Cs(:,:,1)'*(Rs(:,:,1)\Cs(:,:,1));
    F_xw(:,:,1) = [1e4 0; 0 1e12];
    for i=1:N % k=0:w-2
        instA = As(:,:,i);
        instQ = Qs(:,:,i);
        F_xw(:,:,i+1) = -inv(instQ)*instA*inv(F_xw(:,:,i)+instA'*inv(instQ)*instA)*instA'*inv(instQ)+inv(instQ)+Cs(:,:,i+1)'*(Rs(:,:,i+1)\Cs(:,:,i+1));
    end
end