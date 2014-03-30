function X = optimizeSWMmseBcd(K, Q, M, I, J, D, V, L, P)
    X = V;
    for k = 1 : K
        count  = 0;
        prev = 0;
        while true
            count = count + 1;
            if count > 100
                break;
            end
            T = zeros(Q, 1);
            for q = 1 : Q
                rowOffset = (k - 1) * Q;
                colOffset = (k - 1) * I;
                lambdas = L(rowOffset + q, colOffset + 1 : colOffset + I);
                [C, A] = cvector(Q, M, I, D, J, X, lambdas, k, q);
                for i = 1 : I
                    if A(i) == 0
                        rowOffset = (k - 1) * Q * M + (q - 1) * M;
                        colOffset = (k - 1) * I + i;
                        X(rowOffset + 1 : rowOffset + M, colOffset) = zeros(M, 1);
                    end
                end
                if nnz(A) == 0
                    T(q) = 0;
                    continue;
                end
                rowOffset = (k - 1) * Q * M + (q - 1) * M;
                colOffset = (q - 1) * M;
                Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
                rho = spectralRadius(Jkq);
                miuLow = 0;
                miuHigh = upperBoundOfMiu(I, P, A, C);
                miu = (miuLow + miuHigh) / 2;
                delta = zeros(I, 1);
                while true
                    miu = (miuLow + miuHigh) / 2;
                    deltaLow = zeros(I, 1);
                    deltaHigh = upperBoundOfDelta(A, C, I, rho, lambdas, miuHigh);
                    for i = 1 : I
                        if A(i) == 0
                            continue;
                        end
                        t1 = bisectionTarget(Jkq, C(:, i), deltaLow(i), lambdas(i), miu);
                        t2 = bisectionTarget(Jkq, C(:, i), deltaHigh(i), lambdas(i), miu);
                        if (t1 > 1 && t2 > 1) || (t1 < 1 && t2 < 1)
                            fprintf(2, 'parameter selection error\n');
                        end
                        while true
                            delta(i) = (deltaLow(i) + deltaHigh(i)) / 2;
                            target = bisectionTarget(Jkq, C(:, i), delta(i), lambdas(i), miu);
                            if abs(target - 1) < 2e-8
                                break;
                            end
                            if target < 1
                                deltaLow(i) = delta(i);
                            else
                                deltaHigh(i) = delta(i);
                            end
                        end
                    end
                    power = 0;
                    for i = 1 : I
                        if A(i) == 0
                            continue;
                        end
                        power = power + 1 / delta(i)^2;
                    end
                    if abs(miuLow - miuHigh) < 2e-10
                        break;
                    end
                    if power < P
                        miuHigh = miu;
                    else
                        miuLow = miu;
                    end
                end
                for i = 1 : I
                    if A(i) == 0
                        continue;
                    end
                    c = C(:, i);
                    v = (Jkq + (lambdas(i) * delta(i) / 2 + miu) * eye(M)) \ c;
                    rowOffset = (k - 1) * Q * M + (q - 1) * M;
                    colOffset = (k - 1) * I + i;
                    X(rowOffset + 1 : rowOffset + M, colOffset) = v;
                end
                T(q) = miu;
            end
            obj = objectiveSubprolem(Q, M, I, V, J, D, L, k);
            if abs(obj - prev) < 1e-6
                break;
            end
            prev = obj;
        end
    end
    return

function obj = objectiveSubprolem(Q, M, I, V, J, D, L, k)
    obj = 0;
    offset = (k - 1) * Q * M;
    Jk = J(offset + 1 : offset + Q * M, :);
    for i = 1 : I
        rowOffset = (k - 1) * Q * M;
        colOffset = (k - 1) * I + i;
        v = V(rowOffset + 1 : rowOffset + Q * M, colOffset);
        obj = obj + real(v' * Jk * v);
        offset = (k - 1) * I + i;
        d = D(:, offset);
        obj = obj - 2 * real(dot(v, d));
        for q = 1 : Q
            obj = obj + L((k - 1) * Q + q, (k - 1) * I + i) * norm(v((q - 1) * M + 1 : q * M, :));
        end
    end
    return

function y = checkSWMmseConverged(Q, M, I, T, P, V, D, J, L, k)
    y = true;
    for q = 1 : Q
        miu = T(q);
        rowOffset = (k - 1) * Q;
        colOffset = (k - 1) * I;
        lambdas = L(rowOffset + q, colOffset + 1 : colOffset + I);
        [C, A] = cvector(Q, M, I, D, J, V, lambdas, k, q);
        rowOffset = (k - 1) * Q * M + (q - 1) * M;
        colOffset = (q - 1) * M;
        Jkq = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
        for i = 1 : I
            rowOffset = (k - 1) * Q * M + (q - 1) * M;
            colOffset = (k - 1) * I + i;
            v = V(rowOffset + 1 : rowOffset + M, colOffset);
            if norm(v) == 0
                if norm(C(:, i)) <= lambdas(i) / 2
                    continue;
                else
                    y = false;
                    return;
                end
            else
                delta = 1 / norm(v);
                target = bisectionTarget(Jkq, C(:, i), delta, lambdas(i), miu);
                if abs(target - 1) < 2e-13
                    continue;
                else
                    y = false;
                    return;
                end
            end
        end
        power = 0;
        for i = 1 : I
            rowOffset = (k - 1) * Q * M + (q - 1) * M;
            colOffset = (k - 1) * I + i;
            v = V(rowOffset + 1 : rowOffset + M, colOffset);
            power = power + dot(v, v);
        end
        if abs(miu) < 2e-13 || abs(P - power) < 2e-13
            continue;
        else
            y = false;
            return;
        end
    end
    return

function [C, A] = cvector(Q, M, I, D, J, X, lambdas, k, q)
    C = zeros(M, I);
    for i = 1 : I
        rowOffset = (q - 1) * M;
        colOffset = (k - 1) * I + i;
        d = D(rowOffset + 1 : rowOffset + M, colOffset);
        c = d;
        for p = 1 : Q
            if p == q
                continue;
            end
            rowOffset = (k - 1) * Q * M + (q - 1) * M;
            colOffset = (p - 1) * M;
            Jkp = J(rowOffset + 1 : rowOffset + M, colOffset + 1 : colOffset + M);
            rowOffset = (k - 1) * Q * M + (p - 1) * M;
            colOffset = (k - 1) * I + i;
            v = X(rowOffset + 1 : rowOffset + M, colOffset);
            c = c - Jkp * v;
        end
        C(:, i) = c;
    end
    A = zeros(I, 1);
    for i = 1 : I
        if norm(C(:, i), 2) > lambdas(i) / 2
            A(i) = i;
        end
    end
    return

function miuHigh = upperBoundOfMiu(I, P, A, C)
    maxc = 0;
    nz = 0;
    for i = 1 : I
        if A(i) == 0
            continue;
        end
        nz = nz + 1;
        nc = norm(C(:, i));
        if nc > maxc
            maxc = nc;
        end
    end
    miuHigh = (P / nz)^(-0.5) * maxc * 1.1;
    return

function deltaHigh = upperBoundOfDelta(A, C, I, rho, lambdas, miuHigh)
    deltaHigh = zeros(I, 1);
    for i = 1 : I
        if A(i) == 0
            continue;
        end
        deltaHigh(i) = (rho + miuHigh) / (norm(C(:, i)) - lambdas(i) / 2) * 1.1;
    end
    return

function rho = spectralRadius(A)
    rho = max(abs(eig(A)));
    return

function ret = bisectionTarget(Jkq, ciq, delta, lambda, miu)
    ret = delta * norm((Jkq + (lambda * delta / 2 + miu) * eye(size(Jkq))) \ ciq);
    return
