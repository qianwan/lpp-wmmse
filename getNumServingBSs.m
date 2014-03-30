function n = getNumServingBSs(K, Q, M, I, V, reserve)
    n = 0;
    for l = 1 : K
        for q = 1 : Q
            for k = 1 : K
                for i = 1 : I
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    v = V(rowOffset + 1 : rowOffset + M, colOffset);
                    if dot(v, v) > 1.001 * reserve
                        n = n + 1;
                    end
                end
            end
        end
    end
    return
