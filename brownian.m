function [BSs, UEs] = brownian(K, Q, I, locations, outerRadius)
    BSs = zeros(K * Q, 1);
    UEs = zeros(K * I, 1);
    for k = 1 : K
        for q = 1 : Q
            BSs((k - 1) * Q + q) = locations(k);
        end
        for i = 1 : I
            while true
                x = (rand - 0.5) * 2 * outerRadius;
                y = (rand - 0.5) * 2 * outerRadius;
                mx = abs(x);
                my = abs(y);
                valid = true;
                if my > outerRadius * sin(pi / 3) || (sqrt(3) * mx + my > outerRadius * sqrt(3))
                     valid = false;
                end
                if valid == true
                    break;
                end
            end
            UEs((k - 1) * I + i) = x + y * 1j + locations(k);
        end
    end
    return
