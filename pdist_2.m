function d = pdist_2(x,y)
%
%  DESCRIPTION:
%    Pairwise Euclidean distance
%
%  SYNTAX:
%    d = pdist_2(x,y)
%
%  INPUT:
%    x and y are vectors or matrices of same size, where rows are observations 
%    and columns are dimensions.
%
%  OUTPUT:
%    d = distance matrix
%
%  TESTED WITH:
%    Octave 4 and Matlab 2015a
%
%  September, 2018
%  Juan Ignacio Mendoza
%  University of Jyv?skyl?

L = size(x,1);
d = sqrt( repmat(sum(x.^2, 2),1,L) + repmat(sum(y.^2, 2)',L,1) - 2 * x * y') ;

end
