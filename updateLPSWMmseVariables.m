function [U, W, R, obj, Sv] = updateLPSWMmseVariables(K, Q, M, I, N, H, S, V)
    U = zeros(K * I * N, 1);
    W = zeros(K * I, 1);
    R = zeros(K * I, 1);
    Sv = zeros(K * I, 1);
    for ik = 1 : K * I
        rowOffset = (ik - 1) * N + 1 : ik * N;
        hv = H(rowOffset, :) * V;
        C = hv * hv' + eye(N);
        localHv = H(rowOffset, :) * V(:, ik);
        u = C \ localHv;
        U((ik - 1) * N + 1 : ik * N, :) = u;
        W(ik) = 1 / (1 - real(localHv' * u));
        localHvvH = localHv * localHv';
        Lc = C - localHvvH;
        R(ik) = log2(real(det(eye(N) + localHvvH / Lc)));
    end
    obj = sum(R);
    for ik = find(R ~= 0)'
        u = U((ik - 1) * N + 1 : ik * N, 1);
        Sik = S(ik, S(ik, :) ~= 0);
        for ql = Sik
            rowOffset = (ql - 1) * M + 1 : ql * M;
            v = V(rowOffset, ik);
            if norm(v) > 1e-6
                Sv(ik) = Sv(ik) + 1;
            end
        end
    end
    return
