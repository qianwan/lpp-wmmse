function X = optimizeSWMmseFmincon(K, Q, M, I, J, D, V, L, P)
    X = V;
    for k = 1 : K
        offset = (k - 1) * Q * M;
        Jk = J(offset + 1 : offset + Q * M, :);
        offset = (k - 1) * I;
        Dk = D(:, offset + 1 : offset + I);
        rowOffset = (k - 1) * Q * M;
        colOffset = (k - 1) * I;
        lambda = L(k);
        options = optimset('Algorithm', 'interior-point', 'MaxIter', 1500);
        x0 = V(rowOffset + 1 : rowOffset + Q * M, colOffset + 1 : colOffset + I);
        x = fmincon(@(x)fsubproblem(x, Jk, Dk, Q, M, I, lambda), x0, [], [], [], [], [], [], ...
                @(x)powerConstraint(x, Q, M, I, P), options);
        X(rowOffset + 1 : rowOffset + Q * M, colOffset + 1 : colOffset + I) = x;
    end
    return

function obj = fsubproblem(x, Jk, Dk, Q, M, I, lambda)
    obj = 0;
    for i = 1 : I
        obj = obj + real(x(:, i)' * Jk * x(:, i));
        obj = obj - 2 * real(x(:, i)' * Dk(:, i));
        for q = 1 : Q
            obj = obj + lambda * norm(x((q - 1) * M + 1 : q * M, i));
        end
    end
    return

function [c, ceq] = powerConstraint(x, Q, M, I, P)
    ceq = 0;
    c = zeros(Q, 1);
    for q = 1 : Q
        for i = 1 : I
            c(q) = c(q) + x((q - 1) * M + 1 : q * M, i)' * x((q - 1) * M + 1 : q * M, i);
        end
    end
    c = c - P;
    return
