function [U, W, R] = updatePSWMmseVariables(K, Q, M, I, N, H, V)
    U = zeros(K * I * N, 1);
    W = zeros(K * I, 1);
    R = zeros(K * I, 1);
    for k = 1 : K
        for i = 1 : I
            C = zeros(N);
            for l1 = 1 : K
                for l2 = 1 : K
                    rowOffset = (k - 1) * I * N + (i - 1) * N;
                    colOffset = (l1 - 1) * Q * M;
                    h1 = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
                    rowOffset = (l1 - 1) * Q * M;
                    v1 = V(rowOffset + 1 : rowOffset + Q * M, :);
                    rowOffset = (k - 1) * I * N + (i - 1) * N;
                    colOffset = (l2 - 1) * Q * M;
                    h2 = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
                    rowOffset = (l2 - 1) * Q * M;
                    v2 = V(rowOffset + 1 : rowOffset + Q * M, :);
                    C = C + h1 * v1 * v2' * h2';
                end
            end
            C = C + eye(N);
            rowOffset = (k - 1) * I * N + (i - 1) * N;
            h = H(rowOffset + 1 : rowOffset + N, :);
            colOffset = (k - 1) * I + i;
            v = V(:, colOffset);
            localHv = h * v;
            u = C \ localHv;
            rowOffset = (k - 1) * I * N + (i - 1) * N;
            U(rowOffset + 1 : rowOffset + N) = u;
            W((k - 1) * I + i) = 1 / (1 - real(dot(localHv, u)));
            localHvvH = localHv * localHv';
            L = C - localHvvH;
            R((k - 1) * I + i) = log2(real(det(eye(N) + localHvvH / L)));
        end
    end
    return
