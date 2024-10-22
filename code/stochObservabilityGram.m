function [F_x0] = stochObservabilityGram(As,Cs,Qs,Rs,w)
% This function generates the w-step observability Gramian for a
% discrete-time linear time-varying systems with Gaussian noise using the
% recursive formula given in Boyacioglu and van Breugel (2024).

% INPUTS
%   As is the collection n-by-n system matrices
%   Cs is the collection p-by-n output matrices
%   Qs is the collection n-by-n process noise covariance matrices
%   Rs is the collection n-by-n measurement noise covariance matrices
%   w is the number of time steps
% OUTPUTS
%   F_x0 is the stochastic observability Gramian
    F_x0 = Cs(:,:,w)'*(Rs(:,:,w)\Cs(:,:,w));
    for i=1:w-1 % k=0:w-2
        instA = As(:,:,w-i);
        instQ = Qs(:,:,w-i);
        F_x0 = instA'*inv(instQ)*instA-instA'*inv(instQ)*inv(F_x0+inv(instQ))*inv(instQ)*instA+Cs(:,:,w-i)'*(Rs(:,:,w-i)\Cs(:,:,w-i));
    end
end
