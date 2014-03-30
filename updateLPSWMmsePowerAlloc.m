function An = updateLPSWMmsePowerAlloc(K, Q, A, D, S, P, alpha)
    An = A;
    for ql = 1 : K * Q
        Aql = A(:, ql);
        if nnz(Aql) <= 1
            continue;
        end
        Gql = D(:, ql);
        Uql = find(S(:, ql) ~= 0);
        Aql = Aql(Uql);
        Gql = Gql(Uql);
        Aql = Aql - alpha * Gql;
        relaxL = 0;
        relaxH = max(Aql);
        relax = (relaxL + relaxH) / 2;
        proj = waterFilling(relax, Aql);
        while abs(sum(proj) - P) > 1e-4
            if sum(proj) > P
                relaxL = relax;
            elseif sum(proj) < P
                relaxH = relax;
            else
                break;
            end
            relax = (relaxL + relaxH) / 2;
            proj = waterFilling(relax, Aql);
        end
        An(Uql, ql) = proj;
    end
    return

function proj = waterFilling(relax, Aql)
    proj = Aql - relax;
    proj(proj < 0) = 0;
    return
