function L = generateLambdas(K, Q, M, I, N, P, H, SNR, m)
    L = zeros(K * Q, K * I);
    if m == 1
        L = ones(K * Q, K * I) * Q * K / I / sqrt(SNR);
        return;
    elseif m == 2
        for l = 1 : K
            lRowOffset = (l - 1) * Q;
            for q = 1 : Q
                hColOffset = (l - 1) * Q * M + (q - 1) * M;
                for k = 1 : K
                    lColOffset = (k - 1) * I;
                    for i = 1 : I
                        hRowOffset = (k - 1) * I * N + (i - 1) * N;
                        h = H(hRowOffset + 1 : hRowOffset + N, hColOffset + 1 : hColOffset + M);
                        L(lRowOffset + q, lColOffset + i) = model(P, h);
                    end
                end
            end
        end
    end
    return

function ret = model(P, h)
    ret = 1 / sqrt(P * trace(h' * h));
    return
