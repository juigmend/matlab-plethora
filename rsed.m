function d = rsed(a,b)
% 
% DESCRIPTION:   
%   Re-scaled Euclidean distance between two vectors
%
% SYNTAX:
%   d = eud(a,b)
%
% INPUT:
%   a and b are vectors of equal length.
%
% OUTPUT:
%   d: Euclidean Distance
%
% VERSION: 27 June 2020
%
% Juan Ignacio Mendoza
% University of Jyv?skyl?

a = a - min(a);
a = a / max(a);

b = b - min(b);
b = b / max(b);

d = norm(((a-b).^2));

end





