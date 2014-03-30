clear;
K = 4;
M = 4;
N = 2;
Q = 20;
I = 40;
SNRdB = 20;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusters = zeros(K, 1);
r = 2000;
if K == 1
    clusters = 0;
elseif K == 4
    clusters = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
end

numCases = 50;
totalSumRate = 0;
totalNumIterations = 0;
totalNumServingBSs = 0;
maxIterations = inf;
epsilon = 1e-1;
innerEpsilon = 1e-1;
maxNumCand = min(13, floor(K * I / M));
numCand = maxNumCand;
L = ones(K * I, K * Q) * 0;
lassoGamma = 0;
alpha = 3e-1;
innerIter = ones(numCases, 1);
outerIter = zeros(numCases, 1);
rateCase = zeros(numCases, 1);
SvCase = zeros(numCases, 1);

% 0dB, L = 0.015, gamma = 2,  alpha = 0.01,  # = 20, F = 216, S = 0.546
% 0dB, L = 0.02, gamma = 2,   alpha = 0.01,  # = 6, F = 217, S = 0.56
% 0dB, L = 0.05, gamma = 1,   alpha = 0.01,  # = 6, F = 209, S = 0.56
% 0dB, L = 0.05, gamma = 1,   alpha = 0.05,  # = 6, F = 214/212, S = 0.5
% 0dB, L = 0.05, gamma = 1.1, alpha = 0.01,  # = 6, F = 216.8, S = 0.557
% 0dB, L = 0.05, gamma = 1.5, alpha = 0.008, # = 6, F = 212.2/216.2, S = 0.58623
% 0dB, L = 0.05, gamma = 1.5, alpha = 0.01,  # = 6, F = 216.8, S = 0.557
% 0dB, L = 0.05, gamma = 2,   alpha = 0.01,  # = 6, F = 212, S = 0.56
% 0dB, L = 0.05, gamma = 2, alpha = 0.02, # = 6, F = 218.9/211, S = 0.5125/0.52
% 0dB, L = 0.2,  gamma = 1, alpha = 0.01, # = 6, F = 212.8, S = 0.5625
% 0dB, L = 0.2,  gamma = 0.8, alpha = 0.01, # = 6, F = 209.4, S = 0.5688

% 10dB, L = 0.01, gamma = 2,   alpha = 0.01, # = 6, F = 435, S = 4.99
% 10dB, L = 0.01, gamma = 2.5, alpha = 0.01, # = 6, F = 429, S = 4.9
% 10dB, L = 0.01, gamma = 3,   alpha = 0.01, # = 6, F = 440, S = 4.9
% 10dB, L = 0.02, gamma = 1.5, alpha = 0.02, # = 6, F = , S = 

% 0dB,  L = 0.02, gamma = 2,   alpha = 0.01,  # = 6, F = 217, S = 0.56
% 5dB,  L = 0.02, gamma = 2.5, alpha = 0.005, # = 6, F = 318, S = 1.77
% 10dB, L = 0.02, gamma = 2,   alpha = 0.01,  # = 6, F = 435, S = 4.62

% 0dB,  L = 0.02, gamma = 2,   alpha = 0.01,  # = 20, F = 221, S = 0.56

% 0dB,  L = 0, gamma = 0, alpha = 1e-4, # = 50, F = 246.9, S = 5.75
% 0dB,  L = 0, gamma = 0, alpha = 5e-5, # = 50, F = 258.4, S = 5.11

%  0dB, alpha = 1e-3, # = inf, F = 227, S = 0.69, D = 2, C = 1
%  5dB, alpha = 1e-4, # = inf, F = 360, S = 2.7, D = 2, C = 3
%  5dB, alpha = 1e-3, # = inf, F = 361, S = 2.1, D = 2, C = 3
%  5dB, alpha = 1e-2, # = inf, F = 356, S = 1.8, D = 2, C = 3
% 10dB, alpha = 1e-3, # = inf, F = 496, S = 5.1, D = 2, C = 6
% 15dB, alpha = 1e-2, # = 50,  F = 594, S = 6.4, D = 2, C = 11
% 15dB, alpha = 1e-2, # = inf, F = 604, S = 7.8, D = 2, C = 9
% 15dB, alpha = 5e-2, # = inf, F = 607, S = 6.7, D = 2, C = 9
% 20dB, alpha = 1e-1, # = inf, F = 728, S = 10.6, D = 2, C = 12
% 20dB, alpha = 2e-1, # = inf, F = 724, S = 10.1, D = 2, C = 13
% 25dB, alpha = 7e-1, # = inf, F = 767, S = 11.3, D = 2, C = 12
% 30dB, alpha = 5,    # = inf, F = 802, S = 11.7, D = 2, C = 12

for ci = 1 : numCases
    numIterations = 0;
    prev = 0;
    [bss, ues] = brownian(K, Q, I, clusters, r / sqrt(3));
    H = generateMIMOChannel(K, Q, M, bss, I, N, ues, 2);
    numCand = maxNumCand;
    [S, T] = findCandidateBSs(K, Q, M, I, N, H, numCand);
    A = initPowerAllocWithCandidateBSs(K, Q, M, I, H, S, P);
    V = generateRandomTxVectorWithAlloc(K, Q, M, I, A);
    [U, W, R, obj, Sv] = updateLPSWMmseVariables(K, Q, M, I, N, H, S, V, L);
    [mmse, omega, ols] = updateLPSWMmseMatrix(K, Q, M, I, N, H, U, W, S, T, numCand);
    Lp = L ./ (ols.^lassoGamma);
    while abs(prev - obj) > epsilon
        prev = obj;
        numIterations = numIterations + 1;
        if numIterations > maxIterations
            numIterations = numIterations - 1;
            break;
        end
        outerIter(ci) = outerIter(ci) + 1;
        [mmse, omega, ols] = updateLPSWMmseMatrix(K, Q, M, I, N, H, U, W, S, T, numCand);
        % Lp = L ./ (ols.^lassoGamma);
        [V, D, subIter, innerObj, normSt] = updateLPSWMmseTxVector(K, Q, M, I, N, H, A, V, ...
            U, W, S, T, P, Lp, mmse, omega, numCand, 1, false);
        innerPrev = 0;
        innerCnt = 1;
        fprintf(2, '%d.%d.%d obj=%f\n', ci, numIterations, innerCnt, innerObj);
        while abs(innerPrev - innerObj) > innerEpsilon
            innerIter(ci) = innerIter(ci) + 1;
            innerPrev = innerObj;
            An = updateLPSWMmsePowerAlloc(K, Q, A, D, S, P, alpha / innerCnt);
            [X, DT, subIter, innerObj, normSt] = updateLPSWMmseTxVector(K, Q, M, I, N, H, An, V, ...
                U, W, S, T, P, Lp, mmse, omega, numCand, 1, false);
            bAlpha = alpha;
            cc = 1;
            while innerObj > innerPrev
                bAlpha = bAlpha / 2;
                cc = cc + 1;
                if cc > 10
                    % break;
                end
                An = updateLPSWMmsePowerAlloc(K, Q, A, D, S, P, bAlpha / innerCnt);
                [X, DT, subIter, innerObj, normSt] = updateLPSWMmseTxVector(K, Q, M, I, N, H, An, V, ...
                    U, W, S, T, P, Lp, mmse, omega, numCand, 1, false);
            end
            V = X;
            D = DT;
            A = An;
            innerCnt = innerCnt + 1;
            fprintf(2, '%d.%d.%d obj=%f\n', ci, numIterations, innerCnt, innerObj);
        end
        [U, W, R, obj, Sv] = updateLPSWMmseVariables(K, Q, M, I, N, H, S, V, Lp);
        fprintf(2, '%d.%d Obj=%f, R=%f, Sv=%f, SvMin=%d, SvMax=%d\n', ci, numIterations, ...
            obj, sum(R), mean(Sv), min(Sv(Sv ~= 0)), max(Sv));
        [S, T] = updateCandidateBSs(K, Q, M, I, N, H, A, S, numCand);
    end
    rateCase(ci) = sum(R);
    SvCase(ci) = mean(Sv);
    if mod(ci, 10) == 0
        fprintf(2, 'R=%f, Sv=%f\n', mean(rateCase(1 : ci)), mean(SvCase(1 : ci)));
    end
end
fprintf(2, 'R=%f, Sv=%f\n', mean(rateCase), mean(SvCase));
