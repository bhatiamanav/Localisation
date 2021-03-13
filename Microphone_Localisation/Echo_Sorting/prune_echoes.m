function [A, b, V] = prune_echoes(images, loudspeaker, minDistance, As, bs)
%------------------------------------------------------------------------------
% [A, b, V] = prune_echoes(images, loudspeaker, minDistance, As, bs)
%------------------------------------------------------------------------------
%
% Outputs the room vertices based on estimated image sources. Please note
% that this is a simple implementation and likely to have bugs. Stay tuned
% for updated versions.
%
% INPUT :  images           ... estimated image source locations
%          loudspeaker      ... estimated loudspeaker location
%          minDistance      ... minimal distance between image sources and
%                               between vertices
%
% OUTPUT:  A, b             ... plane inequalities describing the room
%                               <A, x> <= b
%          V                ... room vertices
%------------------------------------------------------------------------------


% start with the closest image source

imageDistances = sqrt(sum(bsxfun(@minus, images, loudspeaker).^2));
[~, distanceIndex] = sort(imageDistances);

% keep track of deleted image sources
deleted = zeros(1, numel(imageDistances));

% create the bounding box constraints (arbitrary, con2vert requires this to work)
A = [1  0  0
     0  1  0
     0  0  1
    -1  0  0
     0 -1  0
     0  0 -1; As];
b = [15; 15; 15; 15; 15; 15; bs];

% remove all sources not satisfying the structural input inequalities
if ~isempty(As)
    for s = 1:numel(imageDistances)
        if any(As * images(:, s) > bs)
            deleted(s) = 1;
            fprintf('Discarded IS #: %d (user provided constraint)\n', s);
        end
    end
end

% prune the sources
for i = 1:numel(imageDistances)
    % try expressing this source as a combination of closer sources (for
    % now only second order as a combination of first orders)
    s0 = distanceIndex(i);
    
    for idx1 = 1:(i-1)
        for idx2 = 1:(i-1)
            s1 = distanceIndex(idx1);
            s2 = distanceIndex(idx2);
        


            if (s1 ~= s2) && (~deleted(s1)) && (~deleted(s2))
                % compute an arbitrary point on the "wall" (we don't know if
                % it's a wall yet)
                p2 = (images(:, s2) + loudspeaker) / 2;
                
                % compute the unit, outward pointing normal
                n2 = images(:, s2) - loudspeaker; 
                n2 = n2 / norm(n2);
                
                % combine these two sources
                imageSource12 = images(:, s1) + 2 * ((p2 - images(:, s1))'*n2)*n2;
                
                % do they combine to give the third image source?
                if norm(imageSource12 - images(:, s0)) < minDistance
                    deleted(s0) = 1;
                    fprintf('Discarded IS #: %d (combining)\n', s0);
                end
            end
        end
    end
            
    % if not, see if the corresponding plane interesects the current
    % polyhedron
    
    [V, nr, ~] = lcon2vert(A, b);
    if ~deleted(s0)
        n0 = images(:, s0) - loudspeaker;
        n0 = n0 / norm(n0);
        p0 = (images(:, s0) + loudspeaker) / 2;
        
        A = [A; n0'];
        b = [b; n0'* p0];
        IS(length(b)) = s0;
                
        [V, nrNew, ~] = lcon2vert(A, b);
       
        D = distance(V', V');
        D = D + max(D(:)) * eye(size(D));

        if (numel(nr) == numel(nrNew)) && isempty(setdiff(nrNew, nr) == 0)
            % if the added inequality is redundant, remove it
            
            A = A(1:end-1, :);
            b = b(1:end-1);
            IS = IS(1:end-1);
            fprintf('Discarded IS #: %d (no intersection)\n', s0);
        elseif (min(D(:)) < minDistance / 2)           
            % if some vertices are too close, remove the last plane (heuristic)
            
            A = A(1:end-1, :);
            b = b(1:end-1);
            IS = IS(1:end-1);
            fprintf('Discarded IS #: %d (vertex proximity)\n', s0);
        end
    end
end

[V, nr, ~] = lcon2vert(A, b);
