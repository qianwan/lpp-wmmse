function C = updateCVectors(K, Q, M, I, N, closures, mmse, H, V, U, W)
    C = zeros(K * Q * M, K * I);
    for l = 1 : K
        for q = 1 : Q
            for k = 1 : K
                for i = 1 : I
                    c = zeros(M, 1);
                    for m = closures(l, :)
                        if m == 0
                            continue;
                        end
                        for p = 1 : Q
                            if m == l && p == q
                                continue;
                            end
                            rowOffset = (l - 1) * Q * M + (q - 1) * M;
                            colOffset = (m - 1) * Q * M + (p - 1) * M;
                            mm = mmse(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                            offset = (m - 1) * Q * M + (p - 1) * M;
                            v = V(offset + 1 : offset + M, (k - 1) * I + i);
                            c = c - mm * v;
                        end
                    end
                    w = W((k - 1) * I + i);
                    rowOffset = (k - 1) * I * N + (i - 1) * N;
                    colOffset = (l - 1) * Q * M + (q - 1) * M;
                    h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + M);
                    offset = (k - 1) * I * N + (i - 1) * N;
                    u = U(offset + 1 : offset + N, :);
                    c = c + w * h' * u;
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    C(rowOffset + 1 : rowOffset + M, colOffset) = c;
                end
            end
        end
    end
    return
