function [U, W, R] = updateWMMSEVariables(K, Q, M, I, N, H, V)
    U = zeros(K * I * N, 1);
    W = zeros(K * I, 1);
    R = zeros(K * I, 1);
    for k = 1 : K
        for i = 1 : I
            C = zeros(N);
            for k1 = 1 : K
                rowOffset = (k - 1) * I * N + (i - 1) * N;
                colOffset = (k1 - 1) * Q * M;
                h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
                rowOffset = (k1 - 1) * Q * M;
                v = V(rowOffset + 1 : rowOffset + Q * M, :);
                hv = h * v;
                C = C + hv * hv';
            end
            C = C + eye(N);
            rowOffset = (k - 1) * I * N + (i - 1) * N;
            colOffset = (k - 1) * Q * M;
            h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
            rowOffset = (k - 1) * Q * M;
            colOffset = (k - 1) * I + i;
            v = V(rowOffset + 1 : rowOffset + Q * M, colOffset);
            localHv = h * v;
            offset = (k - 1) * I * N + (i - 1) * N;
            U(offset + 1 : offset + N) = C \ localHv;
            W((k - 1) * I + i) = 1 / (1 - real(dot(U(offset + 1 : offset + N), localHv)));
            L = C - localHv * localHv';
            R((k - 1) * I + i) = log2(real(det(eye(N) + localHv * localHv' / L)));
        end
    end
    return
