function heatwave_duration = DurationHeatwaves(data)
    % Dimensions of the data
    [numLocations,  numDays, numYears] = size(data);
    
    % Initialize the output array
    heatwave_duration = zeros(size(data));
    
    % Loop through each location and year
    for loc = 1:numLocations
        for yr = 1:numYears
            % Get the daily data for this location and year
            dailyData = squeeze(data(loc, :,yr));
            
            % Find sequences of consecutive days exceeding threshold
            consecutiveCount = 0;
            for day = 1:numDays
                if dailyData(day) == 1
                    consecutiveCount = consecutiveCount + 1;
                else
                    if consecutiveCount >= 3
                        heatwave_duration(loc, day-consecutiveCount, yr) = consecutiveCount;
                    end
                    consecutiveCount = 0; % Reset the count
                end
            end
            % Check if the last days of the year form a heatwave
            if consecutiveCount >= 3
                heatwave_duration(loc,end-consecutiveCount+1,yr) = consecutiveCount;
            end
        end
    end
end