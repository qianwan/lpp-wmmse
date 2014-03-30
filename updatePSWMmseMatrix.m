function mmse = updatePSWMmseMatrix(K, Q, M, I, N, H, U, W)
    mmse = zeros(K * Q * M, K * Q * M);
    for k = 1 : K
        for i = 1 : I
            offset = (k - 1) * I * N + (i - 1) * N;
            h = H(offset + 1 : offset + N, :);
            u = U(offset + 1 : offset + N, :);
            w = W((k - 1) * I + i);
            hu = h' * u;
            mmse = mmse + w * (hu * hu');
        end
    end
    return
