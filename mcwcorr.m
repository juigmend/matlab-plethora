function wcorr = mcwcorr(d1,d2,w,h,u)
% Windowed correlation of two one-dimensional signals.
%
% syntax:
%   wcorr = mcwcorr(d1,d2);
%   wcorr = mcwcorr(d1,d2,w);
%   wcorr = mcwcorr(d1,d2,w,h);
%   wcorr = mcwcorr(d1,d2,w,h,u);
%
% input parameters:
%   d1, d2: MoCap data structures, norm strutures or numerical arrays of one column vector
%   w: width of the window (optional, default = 10)
%   h: hop of the window (optional, default = 50)
%   u: width and hop unit, 'percentage' or 'frames' (optional, default = 'percentage')
%
% output:
%   wcorr: a row vector containing the windowed correlation curve
%
% observations:
%   If the length of d1 and d2 are not the same, then zeros will be added at the
%   end of the shortest to match the length of the largest.
%
% VERSION: 6 February 2021
%
% Juan Ignacio Mendoza
% University of Jyväskylä


% process inputs:

if nargin == 1
    disp([10,'There has to be at least two arguments', 10])
    return
end

if nargin == 2
    w = 10;
end

if isstruct(d1)
    data_1 = d1.data;
else
    data_1 = d1;
end

if isstruct(d2)
    data_2 = d2.data;
else
    data_2 = d2;
end

if (size(data_1,2)==1)==0 || (size(data_2,2)==1)==0
    disp([10,'Data has to be uni-dimensional (i.e., one column vector)', 10])
    return
end

length_data_1 = size(data_1,1);
length_data_2 = size(data_2,1);
length_data = length_data_1;

if (length_data_1 == length_data_2) == 0 % check that lengths are equal
    difl = abs(length_data_1 - length_data_2);
    if length_data_1 > length_data_2
        data_2 = vertcat(data_2, zeros(difl,1));
    elseif length_data_1 < length_data_2
        data_1 = vertcat(data_1, zeros(difl,1));
        length_data = length_data_2;
    end
end

if nargin > 2
    if ~isnumeric(w) || ~isnumeric(h)
        disp([10, 'The third and fourth arguments have to be numerical', 10])
        return;
    end
end

if nargin < 5
    u = 'percentage';
    h = 50;
end

if strcmp(u,'percentage')
    window_width = round(length_data * w/100);
    window_hop = round(length_data * h/100);
elseif strcmp(u,'frames')
    window_width = w;
    window_hop = h;
end

margin = floor(window_width / 2);

padding = zeros( margin , 1 );
data_1 = [ padding ; data_1 ; padding ];
data_2 = [ padding ; data_2 ; padding ];

% windowed correlation:
wcorr = zeros(length_data,1);

for i = 1:length_data
    
    this_ending = (window_hop * (i - 1)) + window_width;
    this_starting = this_ending - window_width + 1;
    
    this_window_data_1 = data_1( this_starting : this_ending );
    this_window_data_2 = data_2( this_starting : this_ending );
    
    wcorr(i,1) = corr( this_window_data_1 , this_window_data_2 );  
end

end