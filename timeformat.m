function [HH,MM,SS] = timeformat(s)
%
% syntax:
%   [HH,MM,SS] = timeformat(s)
%
% input:
%   s: seconds
%
% output:
%   HH: hours
%   MM: minutes
%   SS: seconds
%
% VERSION: 8 July 2021
%
% Juan Ignacio Mendoza
% University of Jyväskylä

HH = floor(s / 3600);
s = s - HH * 3600;
MM = floor(s / 60);
SS = s - MM * 60;

end