function d2 = mcorigin(d,rm,DIM)
% Sets coordinates to reference marker at origin, for each frame.
%
% syntax:
%   d2 = mcorigin(d,rm,DIM);
%
% input parameters:
%   d: MoCap structure
%   rm: reference marker
%   DIM: dimensions
%
% output:
%   d2: MoCap structure
%
% example:
%   d2 = mcorigin(d,1,[1,2]);
%   When DIM = [1,2] is the horizontal plane, it Will center each frame horizontally.
%
% VERSION: 9 January 2021
%
% Juan Ignacio Mendoza
% University of Jyväskylä

d2 = d;

ref_cols = [0,0,0];
ref_cols(1:3) = ((rm - 1) * 3 ) + 1;
ref_cols = ref_cols + [0,1,2];

DIM = sort(DIM);

for i_row = 1:size(d.data,1)
    
    for i_dim = DIM
        
        d2.data(i_row,i_dim:3:end) = d2.data(i_row,i_dim:3:end) - d.data(i_row,ref_cols(i_dim));
    end
end

d2.other.origin_marker = rm;