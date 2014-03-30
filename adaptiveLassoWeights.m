function L = adaptiveLassoWeights(K, Q, M, I, V, fun)
    L = zeros(K * Q, K * I);
    for l = 1 : K
        for q = 1 : Q
            for k = 1 : K
                for i = 1 : I
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    v = V(rowOffset + 1 : rowOffset + M, colOffset);
                    L((l - 1) * Q + q, (k - 1) * I + i) = fun(1 / norm(v));
                end
            end
        end
    end
    return
