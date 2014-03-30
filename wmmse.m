clear;
K = 4;
M = 4;
N = 2;
Q = 5;
I = 10;
SNRdB = 0;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
r = 1000;
clusterLocations = [0 + 0j, ...
                    0 + r * 1j, ...
                    r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                    -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
closures = findClusterClosures(clusterLocations, r);
[bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, r / sqrt(3));

numCases = 100;
totalSumRate = 0;
totalNumIterations = 0;
maxIterations = 1e6;
epsilon = 1e-1;
for i = 1 : numCases
    numIterations = 0;
    prev = 0;
    [bsLocations, ueLocations] = brownian(K, Q, I, clusterLocations, r / sqrt(3));
    H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations, 2);
    [V, A] = generateRandomTxVector(K, Q, M, I, N, P, H, closures, 1);
    [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V);
    while abs(prev - sum(R)) > epsilon
        prev = sum(R);
        numIterations = numIterations + 1;
        if numIterations > maxIterations
            numIterations = numIterations - 1;
            break;
        end
        mmse = updateMmseMMatrix(K, Q, M, I, N, H, U, W);
        V = iterateWMMSE(K, Q, M, I, N, mmse, P, H, W, U);
        [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V);
        fprintf(2, 'sum rate @#%d in case#%d: %f\n', numIterations, i , sum(R));
    end
    fprintf(2, 'Case #%d: R = %f, # = %d\n', i, sum(R), numIterations);
    totalSumRate = totalSumRate + sum(R);
    totalNumIterations = totalNumIterations + numIterations;
end
fprintf(2, 'Avg sum rate: %f\n', totalSumRate / numCases);
fprintf(2, 'Avg number of iterations: %f\n', totalNumIterations / numCases);
