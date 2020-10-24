function W = mdu_mask(m, k)
    W = [zeros(m)   ones(m, k) ;
         ones(k, m) zeros(k)  ];
    W = logical(W);
end
%creates the masking matrix for multidimensional unfolding used for
%semidefinite programs
%Given the number of microphones and calibration events it essentially
%returns the mask used for implementing the SNL algorithm.