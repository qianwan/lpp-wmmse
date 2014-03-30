function [S, T] = findCandidateBSs(K, Q, M, I, N, H, maxNumCand)
    S = zeros(K * I, K * Q);
    T = zeros(K * I, maxNumCand * M);
    for ik = 1 : K * I
        TrH = zeros(1, K * Q);
        rowOffset = (ik - 1) * N + 1 : ik * N;
        for ql = 1 : K * Q
            colOffset = (ql - 1) * M + 1 : ql * M;
            h = H(rowOffset, colOffset);
            TrH(ql) = trace(h * h');
        end
        [sTrH, index] = sort(TrH, 'descend');
        S(ik, index(1 : maxNumCand)) = index(1 : maxNumCand);
        Sik = S(ik, S(ik, :) ~= 0);
        for c = 1 : maxNumCand
            T(ik, (c - 1) * M + 1 : c * M) = (Sik(c) - 1) * M + 1 : Sik(c) * M;
        end
    end
    return
