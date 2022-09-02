function [S, d, c, p, l, m] = binseqsi(a,b,L,r)
% 
% DESCRIPTION:   
%   Binary Sequences Similarity version 3
%
% SYNTAX:
%   [S, d, c, p, l, m] = binseqsi(a,b,L,r)
%
% INPUT:
%   a and b are row vectors of equal or different
%   length. These vectors contain indexes of either value in a binary sequence,
%   therefore they should not have repeated elements within each one.
%   L is the length of both binary sequences.
%   r = 0 turns display of results off.
%   r = 1 turns display of results on (default).
%
% OUTPUT:
%       S: Similarity value between 0 and 1;
%          The greater the value, the more similar the sequences are.
%
%       d: Distance of paired elements.
%
%       c: Closeness of paired elements.
%
%       p: Rate of paired elements.
%
%       l: Total lag of paired elements.
%
%       m: Cell array containing the paired elements of a and b.
%
% VERSION: 10 May 2020
%
% Juan Ignacio Mendoza
% University of Jyv?skyl?

% check input arguments:
err_r = ('binseqsi ERROR: r should be 1 or 0');
if exist('r','var')
    if r > 1
        disp(err_r)
        return
    elseif r < 0
        disp(err_r)
        return
    end
else
    r = 1;
end

% check that input vectors do not contain zeroes:
if find([a, b] == 0, 1) > 0
    disp 'binseqsi ERROR: vectors should not contain zeroes.';
    return
end

% check that L is not smaller than the larger of the indexes:
a = sort(a);
b = sort(b);
largerindex = max(a(end),b(end));
if L < largerindex
    disp 'binseqsi ERROR: Length of sequences should not be smaller than the greatest element of both index vectors.';
    return
end
size_a = length(a);
size_b = length(b);

% check that input vectors do not contain duplicates:
errdups = ('binseqsi ERROR: vectors should not contain duplicates.');
if length(unique(a)) ~=  size_a
    disp(errdups);
    return
end
if length(unique(b)) ~=  size_b
    disp(errdups);
    return
end

n_el = size_a + size_b; % total number of elements

% reference matrices for each vector:
id_a = repmat(a',1,size_b);
id_b = repmat(b,size_a,1);

% =========================================================================
% Pair each element of vector b with an element of vector a
% From the paired elements compute Lag, Closeness, Fraction of
% Paired Elements and Similarity.

% initialise matrices:
distmat = zeros(size_a,size_b);
mincols = distmat;
minrows = distmat;

% distance matrix:
for i_1 = 1:size_a
    for i_2 = 1:size_b
        distmat(i_1,i_2) = -(a(i_1)-b(i_2));
    end
end
absdistmat = abs(distmat);

% logical matrix indicating minima of each column:
for i_1 = 1:size_a
    mincols(i_1,:) = absdistmat(i_1,:) == min(absdistmat);
end

% logical matrix indicating minima of each row:
for i_1 = 1:size_b
    minrows(:,i_1) = absdistmat(:,i_1) == min(absdistmat,[],2);
end

allmins = mincols .* minrows; % intersection

n_pel = sum(sum(allmins,1) > 0) + sum(sum(allmins,2) > 0); % total number of paired elements
allmins(allmins == 0) = NaN;

% find to which elements of original a and b do these minima correspond:
from_a = allmins.*id_a;
m{1} =  (from_a(isfinite(from_a)))';
from_b = allmins.*id_b;
m{2} = (from_b(isfinite(from_b)))';

% intermediate measures:
l = nanmean(allmins(:) .* distmat(:)); % lag
d = nansum(allmins(:) .* absdistmat(:)); % distance of paired elements 

% final measures:
c = 1-(d/L); % closeness of paired elements
p = n_pel / n_el; % rate of paired elements
S = c*p;

% =========================================================================
% Display results

if r
        
    disp    --------------------------------------
    disp    '   BINARY SEQUENCES SIMILARITY v3   '
    disp    ' '
    fprintf('                    lag = %f \n', l);
    fprintf('               distance = %f \n', d);
    fprintf('              closeness = %f \n', c);
    fprintf('rate of paired elements = %f \n', p);
    disp ' '
    fprintf('             Similarity = %f \n', S);
    disp ' '
    disp 'Paired elements:'
    disp ' '
    disp(m{1})
    disp(m{2})
    disp    ----------------------------------
    
end

end