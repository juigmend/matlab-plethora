function theta = mctilt(d, m1, m2, m3)
% Returns a vector for the tilt angle of a line from its specified origin.
% This is the angle of the projection of the line on the vertical plane orthogonal
% to the line passing through the origin.
%
% syntax:
%   theta = mctilt(d, m1, m2, m3)
%
% input parameters:
%   d: MoCap structure
%   m1: origin of the line defined by m1 and m2
%   m2: point that defines the angle with the vertical axis
%   m3: define a line passing through m1
%
% output:
%   theta: vector with a tilt angle (radians) for every frame
%
% VERSION: 15 January 2021
%
% Juan Ignacio Mendoza
% University of Jyväskylä

if ~strcmp( d.type , 'MoCap data' ) || ( d.timederOrder ~= 0 )
    
    error('First input shoudl be type = ''MoCap data'' and timederOrder = 0')
end
   
d2 = mcgetmarker(d, [m1,m2,m3]) ;

for i_dim = 1:3
    
    d2.data(:,i_dim:3:end) = d2.data(:,i_dim:3:end) - d2.data(:,i_dim); % translate to origin at m1
end

m3_xy = d2.data( : , 7:8 ); % line on xy
theta_xy = cart2pol( m3_xy(:,1), m3_xy(:,2) ); % azimuth

for i_row = 1:d.nFrames
    
    d2.data(i_row,:) = mcrotate( d2.data(i_row,:) , 180 - ( 180 * theta_xy(i_row) / pi ) , [0 0 1] , [0 0 0] ); % rotate around z
end

m2_yz = d2.data( : , 5:6 ); % line on yz
theta = cart2pol( m2_yz(2), m2_yz(1) ); % tilt