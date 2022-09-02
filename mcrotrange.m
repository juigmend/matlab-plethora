function r = mcrotrange(d, m1, m2, m3, m4)
% Returns a vector for the rotation range of a point m2 around a line described
% by m1 and m3, forming a local coordinate system with the plane defined by m4.
%
% For example, let m1 be the neck, m2 the left shoulder, m3 the root (middle of hips)
% and m4 the left hip. The rotation range will be of the left shoulder respect to
% the neck, orthogonal to the line from the root to the neck.
%
% syntax:
%   r = mcrotrange(d, m1, m2, m3, m4)
%
% input parameters:
%   d: MoCap structure
%   m1: center of rotation
%   m2: rotating point
%   m3 and m4: define a plane for local coordinates
%
% output:
%   r: rotation range (radians)
%
% note:
%   Recommended to use this instead of mcrotationrange.
%
% VERSION: 18 January 2021
%
% Juan Ignacio Mendoza
% University of Jyväskylä

if ~strcmp( d.type , 'MoCap data' ) || ( d.timederOrder ~= 0 )
    
    error('First input shoudl be type = ''MoCap data'' and timederOrder = 0')
end
   
d2 = mcgetmarker(d, [m1,m2,m3,m4]) ;

for i_dim = 1:3
    
    d2.data(:,i_dim:3:end) = d2.data(:,i_dim:3:end) - d2.data(:, i_dim + 6 ); % translate to origin at m3
end

m4_xy = d2.data( : , [10,11] );
theta_xy = cart2pol( m4_xy(:,1), m4_xy(:,2) );

for i_row = 1:d.nFrames
    
    d2.data(i_row,:) = mcrotate( d2.data(i_row,:), -( 180 * theta_xy(i_row) / pi ) , [0 0 1] , [0 0 0] ); % rotate m4 around z
end

m1_xz = d2.data( : , [1,3] );
theta_xz = cart2pol( m1_xz(:,2), m1_xz(:,1) );

for i_row = 1:d.nFrames

    d2.data(i_row,:) = mcrotate( d2.data(i_row,:), -( 180 * theta_xz(i_row) / pi ) , [0 1 0] , [0 0 0] ); % rotate m1 around  y
end

m1_yz_rot = d2.data( :, [2,3] );
phi_yz_rot = cart2pol( m1_yz_rot(:,1) , m1_yz_rot(:,2)  );

for i_row = 1:d.nFrames
    
    d2.data(i_row,:) = mcrotate( d2.data(i_row,:), ( 90 - 180 * phi_yz_rot(i_row) / pi ), [1 0 0] , [0 0 0] ); % rotate m1 around x
end

m2_xy_rot = d2.data( :, [4,5] );
theta = cart2pol( m2_xy_rot(:,1) , m2_xy_rot(:,2)  );
theta_u = unwrap(theta);
r = max(theta_u) - min(theta_u);