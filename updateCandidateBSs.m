function [Sn, T] = updateCandidateBSs(K, Q, M, I, N, H, A, S, numCand)
    Sn = zeros(K * I, K * Q);
    T = zeros(K * I, numCand * M);
    for ik = 1 : K * I
        Aik = A(ik, :);
        Sn(ik, find(Aik ~=0)) = find(Aik ~=0);
        Snik = Sn(ik, :);
        Snik = Snik(Snik ~= 0);
        for c = 1 : nnz(Snik)
            T(ik, (c - 1) * M + 1 : c * M) = (Snik(c) - 1) * M + 1 : Snik(c) * M;
        end

        % Sik = S(ik, :);
        % Sik = Sik(Sik ~= 0);
        % Aik = A(ik, Sik);
        % [sAik, index] = sort(Aik, 'descend');
        % assert(nnz(sAik) <= numCand, 'Alloc should shrink');
        % if nnz(sAik) == numCand
        %     Sn(ik, Sik(index(1 : numCand))) = Sik(index(1 : numCand));
        % elseif nnz(sAik) ~= 0
        %     n = nnz(sAik);
        %     Sn(ik, Sik(index(1 : n))) = Sik(index(1 : n));
        %     rBSs = Sik(index(n + 1 : end));
        %     TrH = zeros(1, length(rBSs));
        %     rowOffset = (ik - 1) * N + 1 : ik * N;
        %     for i = 1 : length(rBSs)
        %         ql = rBSs(i);
        %         colOffset = (ql - 1) * M + 1 : ql * M;
        %         h = H(rowOffset, colOffset);
        %         TrH(i) = trace(h * h');
        %     end
        %     [sTrH, indexTrH] = sort(TrH, 'descend');
        %     Sn(ik, rBSs(indexTrH(1 : numCand - n))) = rBSs(indexTrH(1 : numCand - n));
        % else
        %     continue;
        % end
        % Snik = Sn(ik, :);
        % Snik = Snik(Snik ~= 0);
        % for c = 1 : numCand
        %     T(ik, (c - 1) * M + 1 : c * M) = (Snik(c) - 1) * M + 1 : Snik(c) * M;
        % end
    end
    return
