function closures = findClusterClosures(clusterLocations, radius)
    closures = zeros(length(clusterLocations));
    for k = 1 : size(closures, 1)
    for k1 = 1 : size(closures, 1)
        if abs(clusterLocations(k) - clusterLocations(k1)) <= radius
        closures(k, k1) = k1;
        end
    end
    end
    return
