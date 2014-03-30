function [An, X] = updatePowerAllocation(K, Q, M, I, P, A, S, closures, V, delta, epsilon)
    An = A;
    X = V;
    for l = 1 : K
        for q = 1 : Q
            sub = -S((l - 1) * Q + q, :);
            alloc = A((l - 1) * Q + q, :);
            direc = zeros(1, K * I);
            relaxH = 0;
            for k = closures(l, :)
                if k == 0
                    continue;
                end
                for i = 1 : I
                    a = alloc((k - 1) * I + i);
                    direc((k - 1) * I + i) = a - sub((k - 1) * I + i) * a * delta;
                    if direc((k - 1) * I + i) > relaxH
                        relaxH = direc((k - 1) * I + i);
                    end
                end
            end
            relaxL = 0;
            relaxH = 100;
            relax = (relaxL + relaxH) / 2;
            proj = waterFilling(K, I, closures(l, :), direc, relax, epsilon);
            while abs(sum(proj) - P) > 1e-6
                if sum(proj) > P
                    relaxL = relax;
                elseif sum(proj) < P
                    relaxH = relax;
                else
                    break;
                end
                relax = (relaxL + relaxH) / 2;
                proj = waterFilling(K, I, closures(l, :), direc, relax, epsilon);
            end
            An((l - 1) * Q + q, :) = proj;
            for k = closures(l, :)
                if k == 0
                    continue;
                end
                for i = 1 : I
                    rowOffset = (l - 1) * Q * M + (q - 1) * M;
                    v = V(rowOffset + 1 : rowOffset + M, (k - 1) * I + i);
                    a = An((l - 1) * Q + q, (k - 1) * I + i);
                    if norm(v, 2) == 0
                        v = randn(M, 1) + randn(M, 1) * 1j;
                        X(rowOffset + 1 : rowOffset + M, (k - 1) * I + i) = v / norm(v, 2) * sqrt(a);
                    else
                        X(rowOffset + 1 : rowOffset + M, (k - 1) * I + i) = v / norm(v, 2) * sqrt(a);
                    end
                end
            end
        end
    end
    return

function proj = waterFilling(K, I, closure, direc, relax, reverse)
    proj = zeros(1, K * I);
    for k = closure
        if k == 0
            continue;
        end
        for i = 1 : I
            j = (k - 1) * I + i;
            if direc(j) - relax <= reverse
                proj(j) = reverse;
            else
                proj(j) = direc(j) - relax;
            end
        end
    end
    return
