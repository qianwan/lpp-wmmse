function L = penaltyLPSWMmse(K, Q, M, I, N, H)
    L = zeros(K * I, K * Q);
    for ik = 1 : K * I
        for ql = 1 : K * Q
            rowOffset = (ik - 1) * N + 1 : ik * N;
            colOffset = (ql - 1) * M + 1 : ql * M;
            h = H(rowOffset, colOffset);
            L(ik, ql) = 1 / trace(h * h');
        end
    end
    return
