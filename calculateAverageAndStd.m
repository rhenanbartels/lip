function [averageMass, lowerStd, upperStd] =...
    calculateAverageAndStd(massVector, volumeVector)

    percentageValues = 1:100;
    averageMass = zeros(1, 100);
    lowerStd = zeros(1, 100);
    upperStd = zeros(1, 100);
    
    nPatients = size(massVector, 1);
    
    percentMass = zeros(1, nPatients);
    
    for i = percentageValues
        for j = 1:nPatients
            currentPatientVolume = volumeVector(j, :);
            [val, pos] = min(abs(currentPatientVolume - i));
            percentMass(j) = massVector(j, pos);
        end
        averageMass(i) = mean(percentMass);
        lowerStd(i) = averageMass(i) - 2 * std(percentMass);
        upperStd(i) = averageMass(i) + 2 * std(percentMass);
    end

end