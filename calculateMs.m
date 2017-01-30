function [m3, m15, m85, m97] = calculateMs(massVector, volumeVector)
    nPatients = size(massVector, 1);
    
    m3 = zeros(1, nPatients);
    m15 = zeros(1, nPatients);
    m85 = zeros(1, nPatients);
    m97 = zeros(1, nPatients);
    
    for j = 1:nPatients
        [t, m3Pos] = min(abs(volumeVector(j,:) - 3));
        [t, m15Pos] = min(abs(volumeVector(j,:) - 15));
        [t, m85Pos] = min(abs(volumeVector(j,:) - 85));
        [t, m97Pos] = min(abs(volumeVector(j,:) - 97));
        
        m3(j) = massVector(j, m3Pos);
        m15(j) = massVector(j, m15Pos);
        m85(j) = massVector(j, end) -  massVector(j, m85Pos);
        m97(j) = massVector(j, end) -  massVector(j, m97Pos);
    end
end