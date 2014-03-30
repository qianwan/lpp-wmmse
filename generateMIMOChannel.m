function H = generateMIMOChannel(K, Q, M, bsLocations, I, N, ueLocations, model)
    H = zeros(K * I * N, K * Q * M);
    if model == 1 || nargin == 7
        for i = 1 : size(H, 1)
            for j = 1 : size(H, 2)
                H(i, j) = (randn + randn * 1j) / sqrt(2);
            end
        end
    end
    if model == 2
        for k1 = 1 : K
            for q = 1 : Q
                colOffset = (k1 - 1) * Q * M + (q - 1) * M;
                for k2 = 1 : K
                    for i = 1 : I
                        d = abs(bsLocations((k1 - 1) * Q + q) - ueLocations((k2 - 1) * I + i));
                        rowOffset = (k2 - 1) * I * N + (i - 1) * N;
                        H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + M) = model2(N, M, d);
                    end
                end
            end
        end
    end
    return

function h = model2(N, M, d)
    L = 10^(randn * 8 / 10);
    sigma = sqrt((200 / d)^3 * L);
    h = (randn(N, M) + randn(N, M) * 1j) * sigma;
    return
