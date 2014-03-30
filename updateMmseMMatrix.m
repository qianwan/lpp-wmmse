function mmse = updateMmseMMatrix(K, Q, M, I, N, H, U, W)
    mmse = zeros(K * Q * M, Q * M);
    for k = 1 : K
        m = zeros(Q * M, Q * M);
        for k1 = 1 : K
            for i = 1 : I
                rowOffset = (k1 - 1) * I * N + (i - 1) * N;
                colOffset = (k - 1) * Q * M;
                h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
                offset = (k1 - 1) * I * N + (i - 1) * N;
                u = U(offset + 1 : offset + N, :);
                hu = h' * u;
                m = m + W((k1 - 1) * I + i) * hu * hu';
            end
        end
        offset = (k - 1) * Q * M;
        mmse(offset + 1 : offset + Q * M, :) = m;
    end
    return
