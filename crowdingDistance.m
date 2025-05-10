function distance = crowdingDistance(F, objs)
% F     : Indices of individuals in the current front
% objs  : Objective values matrix (rows: individuals, cols: objectives)
% distance : Crowding distance array for individuals in front F

    numObjs = size(objs, 2);
    numInds = length(F);
    distance = zeros(numInds, 1);
    frontObjs = objs(F, :);  % Objective values of individuals in this front

    for m = 1:numObjs
        % Sort individuals in front by objective m
        [sortedObj, sortIdx] = sort(frontObjs(:, m));
        sortedIdx = F(sortIdx);  % Original indices in population
        
        % Set infinite distance for boundary individuals
        distance(sortIdx(1)) = inf;
        distance(sortIdx(end)) = inf;

        % Normalize objective range
        fMax = max(frontObjs(:, m));
        fMin = min(frontObjs(:, m));
        range = fMax - fMin;

        if range == 0
            continue;  % Skip if all values are the same
        end

        % Compute crowding distance for each individual (except boundaries)
        for i = 2:(numInds - 1)
            d = (sortedObj(i+1) - sortedObj(i-1)) / range;
            distance(sortIdx(i)) = distance(sortIdx(i)) + d;
        end
    end
end