function [S, T] = updateLPSWMmseServBSs(K, Q, M, I, A, maxNumCand)
    S = zeros(K * I, K * Q);
    T = zeros(K * I, maxNumCand * M);
    for ik = 1 : K * I
        Aik = A(ik, :);
        Sik = find(Aik ~= 0);
        S(ik, Sik) = Sik;
        assert(length(Sik) <= maxNumCand, 'new serv BS could not pop up out of nothing');
        for s = 1 : length(Sik)
            T(ik, (s - 1) * M + 1 : s * M) = (Sik(s) - 1) * M + 1 : Sik(s) * M;
        end
    end
    return
