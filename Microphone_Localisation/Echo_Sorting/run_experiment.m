function [estimatedImages, loudspeaker, microphones, estimationErrors] = run_experiment(experimentNumber)
%------------------------------------------------------------------------------
% [estimatedImages, loudspeaker, microphones, estimatedErrors] =
% run_experiment(experimentNumber)
%------------------------------------------------------------------------------
%
% Estimates room shape from impulse responses.
%
% INPUT :  experimentNumber ... number of the experiment (1, 2, 3)
%
% OUTPUT:  estimatedImages  ... estimated coordinates of image sources
%          loudspeaker      ... estimated coordinates of the loudspeaker
%          microphones      ... microphone coordinates (from mic EDM)
%          estimationErrors ... forward error in image source estimation
%------------------------------------------------------------------------------

% get all the parameters of the experiment
[D, directDistances, imageDistances, ~, ~, ~, repeat] = get_experimental_data(experimentNumber);

% TO DO: don't hardcode these
M        = 5;
nPoints  = 22;
thrError = 0.5; % this threshold should be set according to experimental conditions

% setting the last parameter to false permits echo reusing
sortedEchoes = sort_echoes_local(D, imageDistances, 3, max(sqrt(D(:))) * 1.3, nPoints, repeat);

[~, ~, microphones] = dist_opt_mex_3d(D, 3);

nPoints = min(size(sortedEchoes, 1), nPoints);

loudspeaker = trilaterate_beck(microphones, directDistances(:));

estimatedImages  = zeros(3, nPoints);
estimationErrors = zeros(1, nPoints);

fprintf('\n');
fprintf('/----------------------\\\n');
fprintf('|#src | dist/2 | error |\n');
fprintf('|----------------------|\n');
for i = 1:nPoints
    [estimatedImages(:, i) estimationErrors(i)] = trilaterate_beck(microphones, sqrt(sortedEchoes(i, :))');
    fprintf('|%3d  | %5.3f  | %5.3f |\n', i, sqrt(sum(((estimatedImages(:, i) - loudspeaker).^2)))/2, estimationErrors(i));
end
fprintf('\\----------------------/\n\n');

fprintf('Error threshold is set at %f m.\n', thrError);

estimatedImages  = estimatedImages(:, estimationErrors <= thrError);
estimationErrors = estimationErrors(estimationErrors <= thrError);

microphones = microphones';