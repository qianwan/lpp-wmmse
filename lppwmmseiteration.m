clear;
K = 7;
M = 8;
N = 2;
Q = 16;
I = 10;
% 0.125 0.39528 1.25 3.9528 12.5 39.528 125
SNRdB = 20;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusters = zeros(K, 1);
r = 1000;
Dc = r / sqrt(3) * 0.5;
Rc = r / 2 * 1.5;
if K == 1
    clusters = 0;
elseif K == 4
    clusters = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
elseif K == 7
    clusters = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * 1j, ...
                r * cos(pi / 6) + r * (sin(pi / 6) - 1) * 1j, ...
                -r * cos(pi / 6) + r * (sin(pi / 6) - 1) * 1j];
end

numCases = 20;
totalSumRate = 0;
totalNumIterations = 0;
totalNumServingBSs = 0;
maxIterations = 10;
epsilon = 1e-1;
innerEpsilon = 1e-1;
maxNumCand = Q * K;
numCand = maxNumCand;
alpha = 1e-3;
innerIter = ones(numCases, 1);
outerIter = zeros(numCases, 1);

rateCase = zeros(numCases, 1);
SvCase = zeros(numCases, 1);
rpcCase = zeros(numCases, 1);

subIters = [];

for ci = 1 : numCases
    numIterations = 0;
    prev = 0;
    [bss, ues] = brownian(K, Q, I, clusters, r / sqrt(3));
    H = generateMIMOChannel(K, Q, M, bss, I, N, ues, 2);
    numCand = maxNumCand;
    [S, T] = findCandidateBSs(K, Q, M, I, N, H, numCand, bss, ues, Dc, Rc);
    A = initPowerAllocWithCandidateBSs(K, Q, M, I, H, S, P);
    V = generateRandomTxVectorWithAlloc(K, Q, M, I, A);
    [U, W, R, obj, Sv, pc] = updateLPSWMmseVariables(K, Q, M, I, N, H, S, V);
    [mmse, omega] = updateLPSWMmseMatrix(K, Q, M, I, N, H, U, W, S, T, numCand);
    while abs(prev - obj) > epsilon
        prev = obj;
        numIterations = numIterations + 1;
        outerIter(ci) = outerIter(ci) + 1;
        [mmse, omega] = updateLPSWMmseMatrix(K, Q, M, I, N, H, U, W, S, T, numCand);
        [V, D, subIter, innerObj] = updateLPSWMmseTxVector(K, Q, M, I, N, H, A, V, ...
            U, W, S, T, P, mmse, omega, numCand, 12, false);
        subIters = [subIters subIter(subIter~=0)];
        fprintf(2, '#%d\n', ci);
        break;
    end
    rateCase(ci) = sum(R);
    SvCase(ci) = mean(Sv);
    rpcCase(ci) = pc / (P * Q * K);
    if mod(ci, 10) == 0
        fprintf(2, 'R=%f, Sv=%f, rpc=%f\n', mean(rateCase(1 : ci)), mean(SvCase(1 : ci)), mean(rpcCase(1 : ci)));
    end
end
fprintf(2, 'R=%f, Sv=%f, rpc=%f\n', mean(rateCase), mean(SvCase), mean(rpcCase));
