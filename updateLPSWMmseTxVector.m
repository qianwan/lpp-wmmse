function [X, D, Iter, Obj, normSt] = updateLPSWMmseTxVector(K, Q, M, I, N, H, A, V, U, W, S, T, P, Lp, ...
    mmse, omega, maxNumCand, b, standalone)
    X = zeros(size(V));
    D = zeros(K * I, K * Q);
    Iter = zeros(K * I, 1);
    normSt = zeros(K * I, maxNumCand);
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
        subObj = subproblemObjective(Mik, vik, h, w, u, Lp(ik, :), Sik, M);
        prev = 0;
        while abs(prev - subObj) > 1e-4
            Iter(ik) = Iter(ik) + 1;
            prev = subObj;
            % xik = vik; % for parallel
            for index = 1 : length(Sik)
                ql = Sik(index);
                h = H((ik - 1) * N + 1 : ik * N, (ql - 1) * M + 1 : ql * M);
                mmseql = Mik((index - 1) * M + 1 : index * M, :);
                gql = 2 * (mmseql * vik - w * h' * u);
                vql = vik((index - 1) * M + 1 : index * M);
                omegaql = omega(ik, index);
                normvgql = norm(b * omegaql * vql / 2 - gql);
                normSt(ik, index) = normvgql;
                lambda = Lp(ik, ql);
                assert(lambda ~= inf, 'slice error or ols error');
                multiplier = 0;
                nu = A(ik, ql)^(-0.5);
                if normvgql <= lambda
                    vik((index - 1) * M + 1 : index * M) = 0;
                else
                    if A(ik, ql) ~= 0
                        judge = binaryJudge(multiplier, lambda, nu, b, omegaql, normvgql);
                        if judge > 1
                            multiplier = ((normvgql - lambda) * nu - b / 2 * omegaql) / 2;
                        elseif judge < 1
                            nu = b / 2 * omegaql / (normvgql - lambda);
                        end
                        vql = (b * omegaql * vql - 2 * gql) ...
                            / (2 * lambda * nu + 4 * multiplier + b * omegaql);
                        vik((index - 1) * M + 1 : index * M) = vql;
                    else
                        vik((index - 1) * M + 1 : index * M) = 0;
                    end
                end
                D(ik, ql) = -multiplier;
            end
            h = H((ik - 1) * N + 1 : ik * N, Tik);
            subObj = subproblemObjective(Mik, vik, h, w, u, Lp(ik, :), Sik, M);
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

function r = binaryJudge(multiplier, lambda, nu, b, omegaql, normvgql)
    r = nu * 2 * normvgql / (2 * lambda * nu + 4 * multiplier + b * omegaql);
    return

function subObj = subproblemObjective(Mik, v, H, w, u, Lik, S, M)
    subObj = real(v' * Mik * v) - 2 * w * real(u' * H * v);
    for i = 1 : length(S)
        assert(Lik(S(i)) ~= inf, 'penalty should not be inf');
        subObj = subObj + Lik(S(i)) * norm(v((i - 1) * M + 1 : i * M));
    end
    return
