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

numCases = 20;
totalSumRate = 0;
totalNumIterations = 0;
totalNumServingBSs = 0;
maxIterations = 1;
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
subIters = 0;

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
            U, W, S, T, P, Lp, mmse, omega, numCand, 4, false);
        fprintf(2, '%f\n', mean(subIter));
        subIters = subIters + mean(subIter);
        break;
    end
    rateCase(ci) = sum(R);
    SvCase(ci) = mean(Sv);
end
fprintf(2, 'SubIters=%f\n', subIters / numCases);
