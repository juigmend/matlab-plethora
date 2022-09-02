%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                                  FIND MAXIMA                                 %
%                                                                              %
%                                 November 2019                                %
%                                                                              %
%                          Juan Ignacio Mendoza Garay                          %
%                               doctoral student                               %
%                 Department of Music, Art and Culture Studies                 %
%                            University of Jyv?skyl?                           %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program has been tested with:
%   Matlab R2015a
 
% ==============================================================================
% Description:

% Finds indices of maximum values in a vector, if the vector contains one or
% more maximum values.
% The 'max' function alone only finds the index for the first found maximum.

% ==============================================================================
% Instructions:

% If necessary, edit the parameters marked with an arrow like this: <---
% Run the thing, close your eyes and hope for the best.

% ==============================================================================

clc
data = [0, 1, 0, 9, 0, 4, 3, 9, 7, 2]; % <---

max_value = max(data)
max_indexes_1 = find(data == max_value)

% or

max_indexes_2 = find(data ==  max(data))
        