function [X, D, Iter, Obj] = updateLPSWMmseTxVector(K, Q, M, I, N, H, A, V, U, ...
    W, S, T, P, mmse, omega, maxNumCand, b, standalone)
    X = zeros(size(V));
    D = zeros(K * I, K * Q);
    Iter = zeros(1, K * I);
    Obj = 0;
    for ik = 1 : K * I
        u = U((ik - 1) * N + 1 : ik * N);
        if norm(u, 2) == 0
            continue;
        end
        Tik = T(ik, T(ik, :) ~= 0);
        vik = V(Tik, ik);
        Sik = S(ik, S(ik, :) ~= 0);
        offset = (ik - 1) * maxNumCand * M + 1 : (ik - 1) * maxNumCand * M + length(Tik);
        Mik = mmse(1 : length(Tik), offset);
        w = W(ik);
        h = H((ik - 1) * N + 1 : ik * N, Tik);
        subObj = subproblemObjective(Mik, vik, h, w, u, Sik, M);
        prev = 0;
        while abs(prev - subObj) > 1e-6
            prev = subObj;
            if length(Sik) == 48
                Iter(ik) = Iter(ik) + 1;
            end
            xik = vik;
            for index = 1 : length(Sik)
                ql = Sik(index);
                h = H((ik - 1) * N + 1 : ik * N, (ql - 1) * M + 1 : ql * M);
                mmseql = Mik((index - 1) * M + 1 : index * M, :);
                gql = 2 * (mmseql * vik - w * h' * u);
                vql = vik((index - 1) * M + 1 : index * M);
                omegaql = omega(ik, index);
                c = b * omegaql * vql - gql;
                multiplier = 0;
                nu = A(ik, ql)^(-0.5);
                if A(ik, ql) > 1e-7
                    if norm(c, 2) <= b * omegaql * A(ik, ql)^0.5
                        multiplier = 0;
                    else
                        multiplier = (norm(c, 2) * nu - b * omegaql) / 2;
                    end
                    vikql = c / (b * omegaql + 2 * multiplier);
                else
                    vikql = zeros(size(c));
                end
                xik((index - 1) * M + 1 : index * M) = vikql;
                D(ik, ql) = -multiplier;
            end
            vik = xik;
            h = H((ik - 1) * N + 1 : ik * N, Tik);
            subObj = subproblemObjective(Mik, vik, h, w, u, Sik, M);
            if standalone
                assert(subObj <= prev, 'non-monotonic update');
            elseif Iter(ik) > 1
                assert(subObj <= prev, 'non-monotonic update');
            end
        end
        X(Tik, ik) = vik;
        Obj = Obj + subObj;
    end
    return

function subObj = subproblemObjective(Mik, v, H, w, u, S, M)
    subObj = real(v' * Mik * v) - 2 * w * real(u' * H * v);
    return
