function [J, D] = updateSWMmseMatrix(K, Q, M, I, N, H, U, W)
    J = zeros(K * Q * M, Q * M);
    D = zeros(Q * M, K * I);
    for k = 1 : K
        m = zeros(Q * M, Q * M);
        for l = 1 : K
            for j = 1 : I
                w = W((l - 1) * I + j);
                rowOffset = (l - 1) * I * N + (j - 1) * N;
                colOffset = (k - 1) * Q * M;
                h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
                offset = (l - 1) * I * N + (j - 1) * N;
                u = U(offset + 1 : offset + N, :);
                hu = h' * u;
                m = m + w * hu * hu';
            end
        end
        offset = (k - 1) * Q * M;
        J(offset + 1 : offset + Q * M, :) = m;
        for i = 1 : I
            w = W((k - 1) * I + i);
            rowOffset = (k - 1) * I * N + (i - 1) * N;
            colOffset = (k - 1) * Q * M;
            h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
            offset = (k - 1) * I * N + (i - 1) * N;
            u = U(offset + 1 : offset + N, :);
            D(:, (k - 1) * I + i) = w * h' * u;
        end
    end
    return
