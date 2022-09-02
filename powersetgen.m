function P = powersetgen(S, varargin)
% P = powersetgen(S,options)
%
% Version: 4 May 2017
%
% Returns a list P of all subsets of S.
%
% options:
%   'full' includes empty set and S itself
%   'noempty' does not include empty set, but includes S itself (default)
%
% Juan Ignacio Mendoza - 2017
% University of Jyv?skyl?

P{1,1} = NaN;
c = 1;

if nargin >= 2
    if nargin>2 || (xor((strcmp(varargin{1},'full')),strcmp(varargin{1},'noempty')) == 0)
        error('incorrect options')
        return
    elseif strcmp(varargin{1},'full')
        P{1,1} = [];
        c = 2;
    end
end

for i = 1:length(S)
    R{i} = combnk(S,i);
    for i_1 = 1:size(R{i},1)
        P{c,1} = R{i}(i_1,:);
        c = c + 1;
    end
end

end