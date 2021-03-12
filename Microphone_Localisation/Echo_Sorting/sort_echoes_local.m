function echoCombinations = sort_echoes_local(D, echoTimes, dim, windowSizeHalf, nPoints, repeat)
%------------------------------------------------------------------------------
% echoCombinations = sort_echoes(D, echoTimes, dim, windowSizeHalf, nPoints)
%------------------------------------------------------------------------------
%
% Wrapper for 'sort_echoes' that does the local search. That is, it tries
% to combine the peaks from Mic1 only with the peaks that are not too far
% from it. This makes sense for more or less localized microphone arrays.
% The window size parameter should be set to the diameter of the array
% (plus some safety margin to account for the speaker size, etc.).
% 
%
% INPUT : D              ... M by M microphone distance matrix
%         echoTimes      ... M by 1 cell with peak times for each microphone
%         dim            ... ambient dimension (in principle 3)
%         windowSizeHalf ... half of the window size (in meters)
%         nPoints        ... number of highest scoring echo combinations to
%                            return
%
% OUTPUT: echoCombinations ... estimated image source locations for nPoints
%                              highest ranked echo combinations
%         scores           ... scores of the best ranked combinations (lower
%                              is better)
%------------------------------------------------------------------------------


nMicrophones   = numel(echoTimes);
echoTimesLocal = cell(nMicrophones, 1);

echoCombinations  = [];
combinationScores = [];

nEchoesMic1 = length(echoTimes{1});
for i = 1:nEchoesMic1
    echoTimesLocal{1} = echoTimes{1}(i);
    
    noGoodCombinationFlag = false;
    for j = 2:nMicrophones
        echoTimesLocal{j} = echoTimes{j}(abs(echoTimes{j} - echoTimesLocal{1}(1)) <= windowSizeHalf);
        if numel(echoTimesLocal{j}) == 0
            noGoodCombinationFlag = true;
            break;
        end
    end
    if noGoodCombinationFlag
        continue
    end
    
    [combination, score] = sort_echoes(D, echoTimesLocal, dim, 1);
    echoCombinations     = [echoCombinations ; combination];
    combinationScores    = [combinationScores score];
end

[~, bestCombinations] = sort(combinationScores);
nCombinations         = numel(combinationScores);

repeated  = zeros(1, nCombinations);
worstRank = min(nCombinations, 200) + 1;

% For the case when (i-1) doesn't reach nPoints, but it could. For
% instance, when nCombinations = 2.
bestCombinations(worstRank)            = nCombinations + 1; 
echoCombinations(nCombinations + 1, :) = 1e9;

i = worstRank;

%------------------------------------------------------------------------------
% If repeat==true, you allow repeated echoes!

if repeat == false
    for i = 2:worstRank
        if sum(~repeated(1:i-1)) >= nPoints
            break
        end
    
        for j = 1:(i-1)
            if any(echoCombinations(bestCombinations(i), :) == echoCombinations(bestCombinations(j), :))
                repeated(i) = 1;
                break;
            end
        end
    end
end
%------------------------------------------------------------------------------

% Keep only the good ones (non-duplicated)
bestCombinations = bestCombinations(~repeated(1:i-1));
echoCombinations = echoCombinations(bestCombinations, :);


