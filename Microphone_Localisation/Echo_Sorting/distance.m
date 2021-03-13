function d = distance(a,b)

if (nargin ~= 2)
   error('Not enough input arguments');
end

if (size(a,1) ~= size(b,1))
   error('A and B should be of same dimensionality');
end

aa=sum(a.*a,1); bb=sum(b.*b,1); 
d = sqrt(abs(bsxfun(@plus, aa', bb) - 2*a'*b));
