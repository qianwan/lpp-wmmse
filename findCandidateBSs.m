function [S, T] = findCandidateBSs(K, Q, M, I, N, H, maxNumCand, bss, ues, Dc, Rc)
    S = zeros(K * I, K * Q);
    T = zeros(K * I, maxNumCand * M);
    for ik = 1 : K * I
        for ql = 1 : K * Q
            d = abs(bss(ql) - ues(ik));
            if d <= Dc
                S(ik, ql) = ql;
            elseif d < Rc
                S(ik, ql) = ql;
            end
        end
        Sik = S(ik, S(ik, :) ~= 0);
        for c = 1 : maxNumCand
            if c <= length(Sik) && Sik(c) ~= 0
                T(ik, (c - 1) * M + 1 : c * M) = (Sik(c) - 1) * M + 1 : Sik(c) * M;
            end
        end
    end
    return