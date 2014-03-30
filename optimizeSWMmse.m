function V = optimizeSWMmse(K, Q, M, I, J, D, V, L, P, method)
    if strcmp(method, 'cvx')
        V = optimizeSWMmseCvx(K, Q, M, I, J, D, V, L, P);
    elseif strcmp(method, 'bcd')
        V = optimizeSWMmseBcd(K, Q, M, I, J, D, V, L, P);
    elseif strcmp(method, 'fmincon')
        V = optimizeSWMmseFmincon(K, Q, M, I, J, D, V, L, P);
    else
        fprintf(2, '%s: no such method\n', method);
    end
    return
