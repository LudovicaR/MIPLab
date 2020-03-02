function x=demean(x,dim)

% DEMEAN(X) 
% Removes the Average or mean value.
%
% DEMEAN(X,DIM)
% Removes the mean along the dimension DIM of X. 

if(nargin==1),
   dim = 1;
   if(size(x,1) > 1)
      dim = 1;
   elseif(size(x,2) > 1)
      dim = 2;
   end;
end;

dims = size(x);
dimsize = size(x,dim);
dimrep = ones(1,length(dims));
dimrep(dim) = dimsize;

% repmat: replicates the array, according to the number 'dimrep' given, 
% in both dimensions, creating a matrix.
x = x - repmat(mean(x,dim),dimrep);
