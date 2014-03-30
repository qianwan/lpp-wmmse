function [R, servBSs, obj] = deSWMmseUnit(K, Q, M, I, N, H, U, V, W, L, P, numIters, reserve, method)
    R = zeros(K * I, 1);
    obj = 0;
    for i = 1 : numIters
        [J, D] = updateSWMmseMatrix(K, Q, M, I, N, H, U, W);
        V = optimizeSWMmse(K, Q, M, I, J, D, V, L, P, method);
        [U, W, R, obj] = updateSWMmseVariables(K, Q, M, I, N, H, V, L);
    end
    servBSs = getNumServingBSs(K, Q, M, I, V, reserve);
    return
