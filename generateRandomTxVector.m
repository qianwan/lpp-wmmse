function [V, A] = generateRandomTxVector(K, Q, M, I, N, P, H, closures, model)
    V = zeros(K * Q * M, K * I);
    A = zeros(K * Q, K * I);
    for l = 1 : K
        closure = closures(l, :);
        numUEs = nnz(closure) * I;
        power = P / numUEs;
        for q = 1 : Q
            p = zeros(K * I, 1);
            effNoise = ones(K * I, 1) * 1e100;
            for k = closure
                if k == 0
                    continue;
                end
                for i = 1 : I
                    rowOffset = (k - 1) * I * N + (i - 1) * N;
                    colOffset = (l - 1) * Q * M + (q - 1) * M;
                    h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + M);
                    effNoise((k - 1) * I + i) = 1 / trace(h * h');
                end
            end
            p = waterfill(effNoise, P);
            for k = closure
                if k == 0
                    continue;
                end
                for i = 1 : I
                    v = randn(M, 1) + randn(M, 1) * 1j;
                    if model == 1
                        power = power;
                    elseif model == 2
                        power = p((k - 1) * I + i);
                    end
                    v = v / norm(v, 2) * sqrt(power);
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    V(rowOffset + 1 : rowOffset + M, colOffset) = v;
                    A((l - 1) * Q + q, (k - 1) * I + i) = power;
                end
            end
        end
    end
    return
