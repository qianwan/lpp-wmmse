function V = generateRandomTxVectorWithAlloc(K, Q, M, I, A)
    V = zeros(K * Q * M, K * I);
    for ql = 1 : K * Q
        rowOffset = (ql - 1) * M + 1 : ql * M;
        for ik = find(A(:, ql) ~= 0)'
            power = A(ik, ql);
            v = randn(M, 1) + randn(M, 1) * 1j;
            v = v / norm(v, 2) * sqrt(power);
            V(rowOffset, ik) = v;
        end
    end
    return
