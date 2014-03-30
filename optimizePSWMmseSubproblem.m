function [X, S] = optimizePSWMmseSubproblem(K, Q, M, I, N, A, closures, mmse, H, V, U, W, L, reserve)
    X = V;
    S = zeros(K * Q, K * I);
    count = 0;
    maxCount = 70;
    while true
        count = count + 1;
        if count > maxCount
            break;
        end
        fprintf(2, '.');
        C = updateCVectors(K, Q, M, I, N, closures, mmse, H, X, U, W);
        if checkSubproblemConverged(K, Q, M, I, A, C, L, S, X, closures, mmse, reserve)
            break;
        end
        for k = 1 : K
            for i = 1 : I
                for l = closures(k, :)
                    if l == 0
                        continue;
                    end
                    for q = 1 : Q
                        rowOffset = (l - 1) * Q * M + (q - 1) * M;
                        colOffset = (k - 1) * I + i;
                        c = C(rowOffset + 1 : rowOffset + M, colOffset);
                        offset = (l - 1) * Q * M + (q - 1) * M;
                        mm = mmse(offset + 1 : offset + M, offset + 1 : offset + M);
                        a = A((l - 1) * Q + q, (k - 1) * I + i);
                        lambda = L((k - 1) * I + i);
                        [x, multiplier] = blockCoordinateDescent(c, mm, a, lambda, reserve);
                        rowOffset = (l - 1) * Q * M + (q - 1) * M;
                        colOffset = (k - 1) * I + i;
                        X(rowOffset + 1 : rowOffset + M, colOffset) = x;
                        S((l - 1) * Q + q, (k - 1) * I + i) = multiplier;
                    end
                end
            end
        end
    end
    fprintf(2, '\n');
    return
