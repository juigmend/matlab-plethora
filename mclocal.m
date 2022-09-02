function d2 = mclocal(d,rm)
% Sets coordinates to reference marker at origin.
%
% syntax
% d2 = mclocal(d);
%
% input parameters
% d: MoCap structure
% rm: reference marker
%
% output
% d2: MoCap structure
%
%
% VERSION: 9 July 2020
%
% Juan Ignacio Mendoza
% University of Jyv?skyl?

d2 = d;

ref_cols = [0,0,0];
ref_cols(1:3) = ((rm - 1) * 3 ) + 1;
ref_cols = ref_cols + [0,1,2];

for i_row = 1:size(d.data,1)
    
    for i_dim = 1:3
        
        d2.data(i_row,i_dim:3:end) = d2.data(i_row,i_dim:3:end) - d.data(i_row,ref_cols(i_dim));
    end
end

d2.other.origin_marker = rm;