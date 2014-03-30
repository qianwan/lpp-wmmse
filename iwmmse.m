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
clusters = zeros(K, 1);
if K == 1
    clusters = 0 + 0j;
elseif K == 4
    clusters = [0, ...
                r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
end
closures = findClusterClosures(clusters, 0.9 * r);

numCases = 50;
totalSumRate = 0;
totalNumIterations = 0;
maxIterations = 1e6;
epsilon = 1e-1;
TrH = zeros(numCases * K * I * Q, 1);
TrV = zeros(numCases * K * I * Q, 1);
for i = 1 : numCases
    numIterations = 0;
    prev = 0.0;
    [bss, ues] = brownian(K, Q, I, clusters, r / sqrt(3));
    H = generateMIMOChannel(K, Q, M, bss, I, N, ues, 2);
    V = generateRandomTxVector(K, Q, M, I, N, P, H, closures, 1);
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
        fprintf(2, '  %d.%d Sum rate %f\n', i, numIterations, sum(R));
    end
    fprintf(2, '->Case #%d: R = %f, # = %d\n', i, sum(R), numIterations);
    totalSumRate = totalSumRate + sum(R);
    totalNumIterations = totalNumIterations + numIterations;
    fprintf(2, '=>Current avg sum rate: %f\n', totalSumRate / i);
    fprintf(2, '=>Current avg number of iterations: %f\n', totalNumIterations / i);
    for k = 1 : K
        for j = 1 : I
            for q = 1 : Q
                rowOffset = (k - 1) * I * N + (j - 1) * N;
                colOffset = (k - 1) * Q * M + (q - 1) * M;
                h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + M);
                TrH((i - 1) * K * I * Q + (k - 1) * I * Q + (j - 1) * Q + q) = trace(h * h');
                rowOffset = (k - 1) * Q * M + (q - 1) * M;
                colOffset = (k - 1) * I + j;
                v = V(rowOffset + 1 : rowOffset + M, colOffset);
                TrV((i - 1) * K * I * Q + (k - 1) * I * Q + (j - 1) * Q + q) = norm(v);
            end
        end
    end
end
fprintf(2, 'Avg sum rate: %f\n', totalSumRate / numCases);
fprintf(2, 'Avg number of iterations: %f\n', totalNumIterations / numCases);
