function X = optimizeSWMmseCvx(K, Q, M, I, J, D, V, L, P)
    X = zeros(size(V));
    for k = 1 : K
        offset = (k - 1) * Q * M;
        Jk = J(offset + 1 : offset + Q * M, :);
        offset = (k - 1) * I;
        d = D(:, offset + 1 : offset + I);
        lambda = L(k);
        cvx_solver mosek;
        cvx_begin quiet
            cvx_precision best
            variable x(Q * M, I) complex;
            expression part1(I);
            part1(1) = quad_form(x(:, 1), Jk) - 2 * real(x(:, 1)' * d(:, 1));
            for i = 2 : I
                part1(i) = part1(i - 1) + quad_form(x(:, i), Jk) - 2 * real(x(:, i)' * d(:, i));
            end
            expression penalty(I * Q);
            penalty(1) = lambda * norm(x(1 : M, 1));
            for i = 1 : I
                for q = 1 : Q
                    if i == 1 && q == 1
                        continue;
                    end
                    penalty((i - 1) * Q + q) = penalty((i - 1) * Q + q - 1) ...
                        + lambda * norm(x((q - 1) * M + 1 : q * M, i));
                end
            end
            expression mypower(Q, I);
            for q = 1 : Q
                mypower(q, 1) = x((q - 1) * M + 1 : q * M, 1)' * x((q - 1) * M + 1 : q * M, 1);
            end
            for q = 1 : Q
                for i = 2 : I
                    mypower(q, i) = mypower(q, i - 1) + x((q - 1) * M + 1 : q * M, i)' * x((q - 1) * M + 1 : q * M, i);
                end
            end
            minimize(part1(I) + penalty(I * Q));
            subject to
                for q = 1 : Q
                    mypower(q, I) <= P;
                end
        cvx_end
        rowOffset = (k - 1) * Q * M;
        colOffset = (k - 1) * I;
        X(rowOffset + 1 : rowOffset + Q * M, colOffset + 1 : colOffset + I) = x;
    end
    return
