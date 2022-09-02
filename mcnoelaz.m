function d2 = mcnoelaz(d, m1, m2, method)
% Rotates each frame of MoCap data so that a pair of markers remain horizontal (no elevation).
% Optionally, the pair of markers can be locked to remain parallel to the x axis (no azimuth).
%
% syntax:
%   d2 = mc2frontal(d, m1, m2, method);
%
% input parameters:
%   d: MoCap structure
%   m1, m2: numbers of the markers that define the horizontal plane
%   method: 'azlock' locks azimuth to x axis
%           [] will preserve original azimuth
%
% output:
%   d2: MoCap structure
%
% VERSION: 9 January 2021
%
% Juan Ignacio Mendoza
% University of Jyväskylä

if isempty(method)
    reverse_z_rot = true;
elseif strcmp(method,'azlock')
    reverse_z_rot = false;
else
    error('incorrect method')
end

d2 = d;

for frame = 1:size(d.data,1)
    
    tmp = d.data(frame, [ ( 3 * m1 + (-2:0)) , ( 3 * m2 + (-2:0) ) ] );
    
    x1 = tmp(:,1);
    y1 = tmp(:,2);
    z1 = tmp(:,3);
    x2 = tmp(:,4);
    y2 = tmp(:,5);
    z2 = tmp(:,6);
    
    % rotate parallel to z:
    theta = cart2pol( x1-x2 , y1-y2 );
    rot_z = 180-180*theta/pi;
    point_thru_z = [ (x1+x2)/2 , (y1+y2)/2 , 0 ];
    d2.data(frame,:) = mcrotate( d2.data(frame,:), rot_z , [0 0 1], point_thru_z ); 
    
    % rotate parallel to y:
    [~,phi] = cart2sph( x1-x2, y1-y2  , z1-z2 );
    d2.data(frame,:) = mcrotate( d2.data(frame,:), -180*phi/pi, [0 1 0], [ (x1+x2)/2 , 0 , (z1+z2)/2 ] );
    
    % rotate inversely parallel to z:
    if reverse_z_rot
        d2.data(frame,:) = mcrotate( d2.data(frame,:), -rot_z , [0 0 1], point_thru_z );
    end
end
