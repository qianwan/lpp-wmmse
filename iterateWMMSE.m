function V = iterateWMMSE(K, Q, M, I, N, mmse, P, H, W, U)
    V = zeros(K * Q * M, K * I);
    for k = 1 : K
        offset = (k - 1) * Q * M;
        MMatrix = mmse(offset + 1 : offset + Q * M, :);
        power = 0;
        tmp = zeros(Q * M);
        for i = 1 : I
            rowOffset = (k - 1) * I * N + (i - 1) * N;
            colOffset = (k - 1) * Q * M;
            h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
            offset = (k - 1) * I * N + (i - 1) * N;
            u = U(offset + 1 : offset + N, :);
            w = W((k - 1) * I + i);
            v = MMatrix \ (h' * u * w);
            rowOffset = (k - 1) * Q * M;
            colOffset = (k - 1) * I + i;
            V(rowOffset + 1 : rowOffset + Q * M, colOffset) = v;
            power = power + norm(v, 2)^2;
            tmp = tmp + h' * (u * u') * h * w * w;
        end
        if power <= P * Q
            continue
        end
        [D, lambda] = eig(MMatrix);
        phi = D' * tmp * D;
        multiplier = 0;
        miuLow = 0;
        miuHigh = 1;
        while mmseBisectionTarget(phi, lambda, miuHigh) > P * Q
            miuHigh = miuHigh * 2;
        end
        targetValue = 0;
        multiplier = (miuLow + miuHigh) / 2;
        targetValue = mmseBisectionTarget(phi, lambda, multiplier);
        while abs(targetValue - P * Q) / (P * Q) >= 1e-14
            targetValue = mmseBisectionTarget(phi, lambda, multiplier);
            if targetValue > P * Q
                miuLow = multiplier;
            elseif targetValue < P * Q
                miuHigh = multiplier;
            else
                break;
            end
            multiplier = (miuLow + miuHigh) / 2;
        end
        offset = (k - 1) * Q * M;
        MMatrix = mmse(offset + 1 : offset + Q * M, :);
        for i = 1 : I
            rowOffset = (k - 1) * I * N + (i - 1) * N;
            colOffset = (k - 1) * Q * M;
            h = H(rowOffset + 1 : rowOffset + N, colOffset + 1 : colOffset + Q * M);
            offset = (k - 1) * I * N + (i - 1) * N;
            u = U(offset + 1 : offset + N);
            w = W((k - 1) * I + i);
            v = (multiplier * eye(Q * M) + MMatrix) \ h' * u * w;
            rowOffset = (k - 1) * Q * M;
            colOffset = (k - 1) * I + i;
            V(rowOffset + 1 : rowOffset + Q * M, colOffset) = v;
        end
    end
    return

function p = mmseBisectionTarget(phi, lambda, multiplier)
    p = 0;
    for i = 1 : size(phi, 1)
        p = p + real(phi(i, i)) / (real(lambda(i, i)) + multiplier)^2;
    end
return
