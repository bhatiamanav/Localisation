function [echoCombinations, scores] = sort_echoes(D, echoTimes, dim, nPoints)
%---------------------------------------------------------------------------
% [echoCombinations, echoScores] = sort_echoes(D, echoTimes, dim, nPoints)
%---------------------------------------------------------------------------
%
% Approach to "Can One Hear the Shape of a Room" with augmenting the
% Euclidean Distance Matrix with first order echoes.
%
% INPUT : D         ...   M by M microphone distance matrix
%         echoTimes ...   M by 1 cell with peak times for each microphone
%         dim       ...   ambient dimension (in principle 3)
%         nPoints   ...   number of highest scoring echo combinations to
%                         return
%
% OUTPUT: echoCombinations ... estimated image source locations for nPoints
%                              highest ranked echo combinations
%         scores           ... scores of the best ranked combinations (lower
%                              is better)
%---------------------------------------------------------------------------


nMicrophones = numel(echoTimes);

%---------------------------------------------------------------------------
% CONSTRUCT ALL COMBINATIONS OF INPUT ECHOES
%---------------------------------------------------------------------------

lhSide      = '';
rhSide      = '';
rhSideGrid  = '';
for i = 1:nMicrophones
    lhSide      = sprintf('%s t%d'                 ,   lhSide, i);
    rhSide      = sprintf('%s echoTimes{%d}(:).^2,',   rhSide, i);
    rhSideGrid  = sprintf('%s t%d(:) '             ,   rhSideGrid, i);
end

lhSide = ['[' lhSide ']'];
rhSide = ['ndgrid(' rhSide(1:end-1) ')'];

eval([lhSide '=' rhSide ';']);
echoCombinations = []; % to pacify mlint
eval(['echoCombinations = [' rhSideGrid '];']);

%---------------------------------------------------------------------------
% SCORE THE COMBINATIONS ACCORDING TO THE CHOSEN SCORING FUNCTION         
%---------------------------------------------------------------------------

nCombinations     = size(echoCombinations, 1);
combinationScores = zeros(nCombinations, 1);

for i = 1:nCombinations
    if mod(i, 100) == 0
        fprintf('Scoring combination %d/%d\n', i, nCombinations);
    end
    augmentedD           = [D                      echoCombinations(i, :)'; ...
                            echoCombinations(i, :)                      0];
    combinationScores(i) = dist_opt_mex_3d(augmentedD, dim);
end

[~, bestCombinations]    = sort(combinationScores);

%---------------------------------------------------------------------------
% EXCLUDE REUSED ECHOES         
%---------------------------------------------------------------------------

repeated  = zeros(1, nCombinations);
worstRank = min(nCombinations, 200) + 1;

% For the case when (i-1) doesn't reach nPoints, but it could. For
% instance, when nCombinations = 2.
bestCombinations(worstRank)            = nCombinations + 1; 
echoCombinations(nCombinations + 1, :) = 1e9;

for i = 2:worstRank
    if sum(~repeated(1:i-1)) >= nPoints
        break
    end
    
    for j = 1:(i-1)
        if any(echoCombinations(bestCombinations(i), :) == echoCombinations(bestCombinations(j), :))
            repeated(i) = 1;
        end
    end
end

% Keep only the good ones (non-duplicated)
bestCombinations = bestCombinations(~repeated(1:i-1));
echoCombinations = echoCombinations(bestCombinations, :);
scores           = combinationScores(bestCombinations);
