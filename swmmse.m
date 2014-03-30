clear;
K = 4;
M = 4;
N = 2;
Q = 5;
I = 10;
SNRdB = 0;
SNR = 10^(SNRdB / 10);
P = SNR / Q;
clusters = zeros(K, 1);
r = 1000;
if K == 1
    clusters = 0 + 0j;
elseif K == 4
    clusters = [0 + 0j, ...
                0 + r * 1j, ...
                r * cos(pi / 6) + r * sin(pi / 6) * 1j, ...
                -r * cos(pi / 6) + r * sin(pi / 6) * 1j];
end
closures = findClusterClosures(clusters, r * 0.9);
%L = ones(K, 1) * Q * K / I / sqrt(SNR);
[bss, ues] = brownian(K, Q, I, clusters, r / sqrt(3));
H = generateMIMOChannel(K, Q, M, bss, I, N, ues, 2);
L = generateLambdas(K, Q, M, I, N, P, H, SNR, 2);
L = ones(size(L)) * 0.5;
method = 'bcd';

numCases = 50;
totalSumRate = 0;
totalNumIterations = 0;
totalNumServingBSs = 0;
maxIterations = 10;
epsilon = 1e-1;
reserve = 1e-6;
for i = 1 : numCases
    numIterations = 0;
    prev = 0;
    [bss, ues] = brownian(K, Q, I, clusters, r / sqrt(3));
    H = generateMIMOChannel(K, Q, M, bss, I, N, ues, 2);
    [V, A] = generateRandomTxVector(K, Q, M, I, N, P, H, closures, 1);
    [U, W, R, obj] = updateSWMmseVariables(K, Q, M, I, N, H, V, L);
    L = adaptiveLassoWeights(K, Q, M, I, V, @(x)(1.5 * tanh(x / 10)));
    numServgingBSs = 0;
    while abs(prev - obj) > epsilon
        prev = obj;
        numIterations = numIterations + 1;
        if numIterations > maxIterations
            numIterations = numIterations - 1;
            break;
        end
        [J, D] = updateSWMmseMatrix(K, Q, M, I, N, H, U, W);
        V = optimizeSWMmse(K, Q, M, I, J, D, V, L, P, method);
        [U, W, R, obj] = updateSWMmseVariables(K, Q, M, I, N, H, V, L);
        L = adaptiveLassoWeights(K, Q, M, I, V, @(x)(1.5 * tanh(x / 10)));
        numServgingBSs = getNumServingBSs(K, Q, M, I, V, reserve);
        fprintf(2, '  %d.%d Sum rate %f, obj %f, serv BSs %f\n', ...
            i, numIterations, sum(R), obj, numServgingBSs / I / K);
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
