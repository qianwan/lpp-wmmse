function A = initPowerAllocWithCandidateBSs(K, Q, M, I, H, S, P)
    A = zeros(K * I, K * Q);
    for ql = 1 : K * Q
        candUsers = S(:, ql);
        if nnz(candUsers) ~= 0
            powerAlloc = P / nnz(candUsers) * (candUsers ~= 0);
            A(:, ql) = powerAlloc;
        end
    end
    return
