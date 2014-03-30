clear;
K = 4;
M = 4;
N = 2;
Q = 5;
I = 10;
SNRdB = 5;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusterLocations = zeros(1, K);
r = 1000;
if K == 1
    clusterLocations = 0 + 0j;
elseif K == 4
    clusterLocations = [0 + 0j, ...
                        0 + r * 1j, ...
                        r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                        -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
end
closures = findClusterClosures(clusterLocations, r * 1.1);
L = ones(K * I, 1) * 0.5;
[bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, r / sqrt(3));

numCases = 50;
totalSumRate = 0;
totalNumIterations = 0;
totalNumServingBSs = 0;
maxIterations = 1e6;
epsilon = 1e-1;
reserve = 1e-7;
for i = 1 : numCases
    numIterations = 0;
    prev = 0;
    [bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, r / sqrt(3));
    H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations, 2);
    [V, A] = generateRandomTxVector(K, Q, M, I, N, P, H, closures, 1);
    [U, W, rR] = updatePSWMmseVariables(K, Q, M, I, N, H, V);
    R = rR;
    numServgingBSs = 0;
    while abs(sum(prev - rR)) > epsilon
        prev = rR;
        numIterations = numIterations + 1;
        if numIterations > maxIterations
            numIterations = numIterations - 1;
            R = rR;
            break;
        end
        mmse = updatePSWMmseMatrix(K, Q, M, I, N, H, U, W);
        [V, S] = optimizePSWMmseSubproblem(K, Q, M, I, N, A, closures, mmse, H, V, U, W, L, reserve);
        [U, W, rR] = updatePSWMmseVariables(K, Q, M, I, N, H, V);
        numServgingBSs = getNumServingBSs(K, Q, M, I, V, reserve);
        fprintf(2, '  %d.%d Sum rate %f, serv BSs %f\n', i, numIterations, sum(rR), numServgingBSs / I / K);
        if sum(rR - prev) < epsilon
            R = prev;
            numIterations = numIterations - 1;
            break;
        end
        [A, V] = updatePowerAllocation(K, Q, M, I, P, A, S, closures, V, (sum(rR - prev)) * 0.025, reserve);
        [U, W, R] = updatePSWMmseVariables(K, Q, M, I, N, H, V);
    end
    totalSumRate = totalSumRate + sum(R);
    totalNumIterations = totalNumIterations + numIterations;
    totalNumServingBSs = totalNumServingBSs + numServgingBSs;
    fprintf(2, '->Case #%d: R = %f # = %d\n', i, sum(R), numIterations);
    fprintf(2, '=>Current avg sum rate: %f\n', totalSumRate / i);
    fprintf(2, '=>Current avg number of iterations: %f\n', totalNumIterations / i);
    fprintf(2, '=>Current avg number of serving BSs per user: %f\n', totalNumServingBSs / i / K / I);
end
fprintf(2, 'Avg sum rate: %f\n', totalSumRate / numCases);
fprintf(2, 'Avg number of iterations: %f\n', totalNumIterations / numCases);
fprintf(2, 'Avg number of serving BSs per user: %f\n', totalNumServingBSs / numCases / K / I);
