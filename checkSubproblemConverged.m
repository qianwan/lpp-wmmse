function converged = checkSubproblemConverged(K, Q, M, I, A, C, L, S, X, closures, mmse, reserve)
    converged = true;
    for l = 1 : K
        for q = 1 : Q
            for k = closures(l, :)
                if k == 0
                    continue;
                end
                for i = 1 : I
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    c = C(rowOffset + 1 : rowOffset + M, colOffset);
                    lambda = L((k - 1) * I + i);
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    x = X(rowOffset + 1 : rowOffset + M, colOffset);
                    a = A((l - 1) * Q + q, (k - 1) * I + i);
                    if norm(c, 2) <= lambda / 2 || a <= 1.001 * reserve
                        if norm(x, 2) ~= 0
                            converged = false;
                            return
                        end
                    elseif norm(x) > 0
                        multiplier = S((l - 1) * Q + q, (k - 1) * I + i);
                        offset = (l - 1) * Q * M + (q - 1) * M;
                        mm = mmse(offset + 1 : offset + M, offset + 1 : offset + M);
                        if abs(multiplier * (a - norm(x, 2)^2)) > 1e-3
                            converged = false;
                            return
                        end
                        rep = -2 * norm(x, 2) / lambda * (multiplier * x + mm * x - c);
                        if norm(rep - x, 2)^2 > (0.01 * norm(x, 2)^2)
                            converged = false;
                            return
                        end
                    else
                        converged = false;
                        return
                    end
                end
            end
        end
    end
    return
