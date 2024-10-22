% This code generates the figure seen in the paper titled "Duality of
% Observability and Constuctability and Their Relation to the Fisher
% Information" by Burak Boyacioglu and Floris van Breugel. It is written by
% Burak Boyacioglu in October 2024.

N = 30;

A = [2 -1; 0 1];
C = [1 0];

Q = 1e-2*[3.6 1.2; 1.2 6];
R = 0.1;

n = size(A,2); % no. of state variables
p = size(C,1); % no. of outputs

As = zeros(n,n,N+1);
Cs = zeros(p,n,N+1);
Qs = zeros(n,n,N+1);
Rs = zeros(p,p,N+1);

for k=0:N
    As(:,:,k+1) = A+[0 sin(pi*k/18); cos(pi*k/18) 0]; % A is time-varying
    Cs(:,:,k+1) = C;
    Qs(:,:,k+1) = Q;
    Rs(:,:,k+1) = R;
end

ws = 1:N+1; % window sizes

% MAIN METHOD
FIM0_recursive = zeros(n,n,length(ws)); 

FIM0_recursive(:,:,1) =  Cs(:,:,1)'*(Rs(:,:,1)\Cs(:,:,1));
for w=2:N+1   
    FIM0_recursive(:,:,w) = stochObservabilityGram(As,Cs,Qs,Rs,w);
end

% METHOD GIVEN IN TENNY AND RAWLINGS (2002)
FIM0_TnR = zeros(n,n,length(ws));
FIM0_TnR(:,:,1) =  Cs(:,:,1)'*(Rs(:,:,1)\Cs(:,:,1));

for w=2:N+1
    
    O = Cs(:,:,1);
    prodA = eye(n);
    for i = 1:(w-1)
        prodA = As(:,:,i)*prodA;
        O = [O; Cs(:,:,i+1)*prodA]; % Append the product of the last p rows of O and A
    end
    
    M = zeros(w*p, w*n);
    
    for blockRow = 1:w
        for blockCol = 1:w
            rowIndex = (blockRow-1)*p + (1:p); % Calculate row index for block
            colIndex = (blockCol-1)*n + (1:n); % Calculate column index for block
            if blockRow == 1
                % First row of blocks is all zeros
                continue;
            elseif blockRow-blockCol == 1
                M(rowIndex, colIndex) = Cs(:,:,blockRow);
            elseif blockRow-blockCol > 1
                M(rowIndex, colIndex) = Cs(:,:,blockRow);
                for i=blockRow-1:-1:blockCol+1
                    M(rowIndex, colIndex) = M(rowIndex, colIndex)*As(:,:,i);
                end
            end
        end
    end
    q = num2cell(Qs(:,:,1:w),[1,2]);
    blkdiagQ = blkdiag(q{:});
    r = num2cell(Rs(:,:,1:w),[1,2]);
    blkdiagR = blkdiag(r{:});

    bigR = blkdiagR+M*blkdiagQ*M';
    FIM0_TnR(:,:,w) = O'*inv(bigR)*O;
%     FIM0_TnR_noInv(:,:,w) = O'*(bigR\O); % numerically more stable, but not good enough
end

% METHOD GIVEN IN THEOREM 1
FIM0_thm1 = zeros(n,n,length(ws));
FIM0_thm1(:,:,1) =  Cs(:,:,1)'*(Rs(:,:,1)\Cs(:,:,1));

PHIs = zeros(n,n,N,N); % construct state transition matrices
for i = 1:N
    PHIs(:,:,i,i) = eye(n);
    for j = i+1:N
        PHIs(:,:,j,i) = As(:,:,j)*PHIs(:,:,j-1,i); % \Phi_{j,i}
    end
end

for w=2:N+1
    
    O = Cs(:,:,1);
    prodA = eye(n);
    for i = 1:(w-1)
        prodA = As(:,:,i)*prodA;
        O = [O; Cs(:,:,i+1)*prodA]; % Append the product of the last p rows of O and A
    end
    
    subR = zeros(p*(w-1), p*(w-1));
    
    for j = 1:w-1
        for k = 1:w-1
            % Calculate R_{jk}
            if j == k
                subR((j-1)*p+1:j*p, (k-1)*p+1:k*p) = Rs(:,:,j+1);
                for i=1:k
                    subR((j-1)*p+1:j*p, (k-1)*p+1:k*p) = subR((j-1)*p+1:j*p, (k-1)*p+1:k*p)+Cs(:,:,j+1)*PHIs(:,:,j,i)*Qs(:,:,i)*PHIs(:,:,j,i)'*Cs(:,:,j+1)';
                end
            elseif j > k
                for i=1:k
                    subR((j-1)*p+1:j*p, (k-1)*p+1:k*p) = subR((j-1)*p+1:j*p, (k-1)*p+1:k*p)+Cs(:,:,j+1)*PHIs(:,:,j,i)*Qs(:,:,i)*PHIs(:,:,k,i)'*Cs(:,:,k+1)';
                end
                % Assign the corresponding upper triangle element
                subR((k-1)*p+1:k*p, (j-1)*p+1:j*p) = subR((j-1)*p+1:j*p, (k-1)*p+1:k*p)';
            end
        end
    end
    bigR = [Rs(:,:,1) zeros(p,p*(w-1)); zeros(p*(w-1),p) subR];
    FIM0_thm1(:,:,w) = O'*inv(bigR)*O;
%     FIM0_noInv(:,:,w) = O'*(bigR\O); % numerically more stable, but not good enough
end

%% 11- and 33-step Dual Constructability Calculations (for illustrative purposes)
w = 11;
F0_N10 = zeros(n,n,w);
F0_N10(:,:,1) = Cs(:,:,w)'*(Rs(:,:,w)\Cs(:,:,w));
for i=1:w-1 % k=0:w-2
    instA = As(:,:,w-i); instQ = Qs(:,:,w-i);
    F0_N10(:,:,i+1) = instA'*inv(instQ)*instA-instA'*inv(instQ)*inv(F0_N10(:,:,i)+inv(instQ))*inv(instQ)*instA+Cs(:,:,w-i)'*(Rs(:,:,w-i)\Cs(:,:,w-i));
end

w = 31;
F0_N30 = zeros(n,n,w);
F0_N30(:,:,1) = Cs(:,:,w)'*(Rs(:,:,w)\Cs(:,:,w));
for i=1:w-1 % k=0:w-2
    instA = As(:,:,w-i); instQ = Qs(:,:,w-i);
    F0_N30(:,:,i+1) = instA'*inv(instQ)*instA-instA'*inv(instQ)*inv(F0_N30(:,:,i)+inv(instQ))*inv(instQ)*instA+Cs(:,:,w-i)'*(Rs(:,:,w-i)\Cs(:,:,w-i));
end

%% FIGURES

dualityGreen = [106,180,46]/256;

% Create figure with specified size
figure('Units', 'centimeters', 'Position', [2, 2, 18, 18]);

% Set up a 2x2 tiled layout for the four panels, adding padding for the legend
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% First subplot
nexttile;
h1 = plot(ws, squeeze(FIM0_recursive(1,1,:)), 'k', 'LineWidth', 4.0);
hold on;
h2 = plot(ws, squeeze(FIM0_thm1(1,1,1:N+1)), 'LineWidth', 1.5);
h3 = plot(ws, squeeze(FIM0_TnR(1,1,1:N+1)), '--', 'LineWidth', 1.5);
h4 = plot(1:11, squeeze(F0_N10(1,1,:)), ':', 'Color', dualityGreen, 'LineWidth', 1.5);
h5 = plot(ws, squeeze(F0_N30(1,1,:)), 'Color', dualityGreen, 'LineWidth', 1.5);
ylabel('$[F_{\downarrow_{w}}^{\mathbf{x}_{0}}]_{1,1}$', 'Interpreter', 'latex');

% Second subplot
nexttile;
plot(ws, squeeze(FIM0_recursive(1,2,:)), 'k', 'LineWidth', 4.0);
hold on;
plot(ws, squeeze(FIM0_thm1(1,2,1:N+1)), 'LineWidth', 1.5);
plot(ws, squeeze(FIM0_TnR(1,2,1:N+1)), '--', 'LineWidth', 1.5);
plot(1:11, squeeze(F0_N10(1,2,:)), ':', 'Color', dualityGreen, 'LineWidth', 1.5);
plot(ws, squeeze(F0_N30(1,2,:)), 'Color', dualityGreen, 'LineWidth', 1.5);
ylabel('$[F_{\downarrow_{w}}^{\mathbf{x}_{0}}]_{1,2}$', 'Interpreter', 'latex');

% Third subplot
nexttile;
plot(ws, squeeze(FIM0_recursive(2,1,:)), 'k', 'LineWidth', 4.0);
hold on;
plot(ws, squeeze(FIM0_thm1(2,1,1:N+1)), 'LineWidth', 1.5);
plot(ws, squeeze(FIM0_TnR(2,1,1:N+1)), '--', 'LineWidth', 1.5);
plot(1:11, squeeze(F0_N10(2,1,:)), ':', 'Color', dualityGreen, 'LineWidth', 1.5);
plot(ws, squeeze(F0_N30(2,1,:)), 'Color', dualityGreen, 'LineWidth', 1.5);
ylabel('$[F_{\downarrow_{w}}^{\mathbf{x}_{0}}]_{2,1}$', 'Interpreter', 'latex');
xlabel('$\mathrm{window~size},\ w$', 'Interpreter', 'latex');

% Fourth subplot
nexttile;
plot(ws, squeeze(FIM0_recursive(2,2,:)), 'k', 'LineWidth', 4.0);
hold on;
plot(ws, squeeze(FIM0_thm1(2,2,1:N+1)), 'LineWidth', 1.5);
plot(ws, squeeze(FIM0_TnR(2,2,1:N+1)), '--', 'LineWidth', 1.5);
plot(1:11, squeeze(F0_N10(2,2,:)), ':', 'Color', dualityGreen, 'LineWidth', 1.5);
plot(ws, squeeze(F0_N30(2,2,:)), 'Color', dualityGreen, 'LineWidth', 1.5);
ylabel('$[F_{\downarrow_{w}}^{\mathbf{x}_{0}}]_{2,2}$', 'Interpreter', 'latex');
xlabel('$\mathrm{window~size},\ w$', 'Interpreter', 'latex');

% Collect all plot handles
allHandles = [h1, h2, h3, h4, h5];

% Add a legend below all plots
lgd = legend(allHandles, {'recursive formula', 'Theorem 1', 'Tenny and Rawlings (2002)','11-step dual system''s constructability', '31-step dual system''s constructability'}, ...
    'Interpreter', 'latex', 'Location', 'southoutside', 'NumColumns', 3);

% Adjust legend to fit below all plots
lgd.Layout.Tile = 'south';
