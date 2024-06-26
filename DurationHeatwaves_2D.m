function heatwave_duration = DurationHeatwaves_2D(data)
    % Dimensions of the data
    [numLocations, numDays] = size(data);
    
    % Initialize the output array
    heatwave_duration = zeros(size(data));
    
    % Loop through each location
    for loc = 1:numLocations
        % Get the daily data for this location
        dailyData = data(loc, :);
        
        % Find sequences of consecutive days exceeding threshold
        consecutiveCount = 0;
        for day = 1:numDays
            if dailyData(day) == 1
                consecutiveCount = consecutiveCount + 1;
            else
                if consecutiveCount >= 3
                    heatwave_duration(loc, day - consecutiveCount) = consecutiveCount;
                end
                consecutiveCount = 0; % Reset the count
            end
        end
        % Check if the last days in the data form a heatwave
        if consecutiveCount >= 3
            heatwave_duration(loc, end - consecutiveCount + 1) = consecutiveCount;
        end
    end
end
