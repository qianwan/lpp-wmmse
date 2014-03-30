function [x, multiplier] = blockCoordinateDescent(c, mm, a, lambda, reserve)
    M = size(mm, 1);
    if norm(c, 2) <= lambda / 2 || a < 0.9 * reserve
        x = zeros(M, 1);
        multiplier = 0;
    else
        theta = 1 / sqrt(a);
        multiplier = 0;
        target = bisectionTarget(0, theta, lambda, mm, c);
        if target > 1
            miuLow = 0;
            miuHigh = norm(c) * theta * 1.1;
            multiplier = (miuLow + miuHigh) / 2;
            theta = 1 / sqrt(a);
            target = bisectionTarget(multiplier, theta, lambda, mm, c);
            while abs(target - 1) > 1e-6
                if target > 1
                    miuLow = multiplier;
                elseif target < 1
                    miuHigh = multiplier;
                end
                multiplier = (miuLow + miuHigh) / 2;
                target = bisectionTarget(multiplier, theta, lambda, mm, c);
            end
        elseif target < 1
            tLow = 0;
            tHigh = (spectralRadius(mm) + norm(c, 2) * theta * 1.1) / (norm(c, 2) - lambda / 2);
            multiplier = 0;
            theta = (tLow + tHigh) / 2;
            target = bisectionTarget(multiplier, theta, lambda, mm, c);
            while abs(target - 1) > 1e-6
                if target > 1
                    tHigh = theta;
                elseif target < 1
                    tLow = theta;
                end
                theta = (tLow + tHigh) / 2;
                target = bisectionTarget(multiplier, theta, lambda, mm, c);
            end
        end
        x = (mm + (lambda * theta / 2 + multiplier) * eye(size(mm))) \ c;
    end
    return

function rho = spectralRadius(A)
    rho = max(abs(eig(A)));
    return

function r = bisectionTarget(multiplier, theta, lambda, mm, c)
    r = theta * norm((mm + (lambda * theta / 2 + multiplier) * eye(size(mm))) \ c);
    return
