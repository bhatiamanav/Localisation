====================================
| CAN ONE HEAR THE SHAPE OF A ROOM |
====================================

This package contains Matlab and C (mex) files accompanying the paper 
"Acoustic echoes reveal room shape" by I. Dokmanic, R. Parhizkar, A. Walther,
Y. M. Lu and M. Vetterli.

Note that the package contains a mex file dist_opt_mex_3d.c, that should be 
compiled on your platform prior to using the code. For more information on mex
files:

>> help mex




=========
| FILES |
=========

The main echo-sorting routine is sort_echoes_local(). It calls another routine,
sort_echoes(), that tries all the combinations of input echo times, which could
be slow. Therefore sort_echoes_local() only tries "reasonable" combinations, 
taking into account the size of the microphone array.




==============
| HOW TO RUN |
==============

To get the "sorted" echoes, do the following

>> [images, loudspeaker, microphones, ~] = run_experiment(experiment);

experiment = 1: omnidirectional loudspeaker, classroom
experiment = 2: directional speaker, classroom
experiment = 3: cathedral portal


From these sorted echoes we can get the room geometry by discarding the
higher-order image sources. To do so, run the following

>> [A, b, V] = prune_echoes(images, loudspeaker, 2, As, bs);

For Experiment 1, set As = [], bs = []. For Experiment 2 we used a directional 
loudspeaker, and we don't hope to reconstruct the wall directly behind the 
loudspeaker. We can take this into account by setting:

>> As = (loudspeaker - microphones(:, 5))' / norm(loudspeaker - microphones(:, 5)); 
>> bs = dot(As', loudspeaker + As' * 0.1);

This will communicate to prune_echoes() to disregard any image source "x" 
not satisfying the inequality <As, x> <= bs. The inequality is setup so that 
it roughly means "discard everything behind the loudspeaker".

To visualize the result you can use the following code:

>> figure; clf; hold all;
>> scatter3(V(:, 1), V(:, 2), V(:, 3), 'r*'); % room vertices
>> scatter3(microphones(1, :), microphones(2, :), microphones(3, :), 'k>', 'filled');
>> scatter3(loudspeaker(1), loudspeaker(2), loudspeaker(3), 'bo', 'filled');
>> view(3); axis equal; grid on;
>> legend('Room vertices', 'Microphones', 'Loudspeaker');

Note that the result will be correct up to a rotation and a reflection. We 
could rotate and mirror it into an "absolutely correct" orientation if we
assume, for example, the knowledge of the floor reflection and the microphone
array orientation.

In Experiment 3, as described in the paper, the "room" does not satisfy the
assumptions of it being a convex polyhedron, thus it makes no sense to use the 
function prune_echoes(). The following distances correspond to the dimensions
of the portal:

>> norm(images(:, 6) - images(:, 20))/2

ans =

    6.7587

>> norm(loudspeaker - images(:, 22))/2

ans =

    6.2751

To this last distance, we must to add the distance of the acoustic center of 
the loudspeaker to the back wall (a bit more than 20 cm) to get the dimension.
In terms of angles, walls 6 and 20 are esentially parallel, but wall 22 is not
exactly perpendicular to them. A probable reason is that in echo sorting,
echoes from wall 22 received by 4 microphones were combined with a different 
echo off an obstacle received by the 5th microphone (since the true one is 
missing).



====================
| ACKNOWLEDGEMENTS |
====================

This code uses functions lcon2vert and vert2lcon available at:

http://www.mathworks.com/matlabcentral/fileexchange/30892-representing-polyhedral-convex-hulls-by-vertices-or-inequalities

And the function distance available at:

http://www.mathworks.com/matlabcentral/fileexchange/71-distance-m



