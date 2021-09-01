function u = getcartesiandirection(ori,loc);
%  GETCARTESIANDIRECTION -- Construct Cartesian unit vector from
%                           Earth-based orientation
%
%  getcartesiandirection(ori,loc) calculates the unit vector u in
%  Earth-based Cartesian coordinates corresponding to an orientation
%  described by the structure ori, at a location on the Earth
%  described by the structure loc.
%
%  The output is a unit vector u in an Earth-fixed Cartesian
%  coordinate system whose first, second, and third axes pierce the
%  surface of the Earth at 
%         1. The intersection of the Equator and the Prime Meridian
%         2. The intersection of the Equator and the meridian at 90
%            degrees East longitude
%         3. The North Pole
%
%  The input structure loc has the fields
%        lat: geodetic latitude (measured North from the Equator) in radians
%        lon: geodetic longitude (measured East from the Prime
%             Meridian) in radians
%     height: elevation in meters above the WGS-84 reference ellipsoid 
%  Such a structure can be created from the geographic coordinates (in
%  degrees) using the function CREATELOCATION
%
%  The input orientation is a structure ori with the fields
%        az: azimuth in radians East (clockwise) of North
%       alt: altitude (tilt) angle in radians above the local tangent plane
%  Such a structure can be created from the local angles (in degrees)
%  using the function CREATEORIENTATION
%  
%  The function first calculates the components of the unit vector in
%  a local Cartesian basis (with axes aligned in the local East, North,
%  and Up directions), then transforms the vector first into an
%  orthonormal basis based on Earth-fixed cylindrical coordinates and
%  finally into Earth-fixed Cartesian coordinates.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%  $Id: getcartesiandirection.m,v 1.4 2006/04/14 21:05:42 whelan Exp $
%
%  See also CREATELOCATION, CREATEORIENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Components of u in a local Cartesian basis of East, North, Up

uloc = transpose([cos(ori.alt)*sin(ori.az) ...
                  cos(ori.alt)*cos(ori.az) ...
                  sin(ori.alt) ]);

coslat = cos(loc.lat);
sinlat = sin(loc.lat);
coslon = cos(loc.lon);
sinlon = sin(loc.lon);

% Transformation matrix from E,N,U to cylindrical basis
Rcylloc = [ 0 -sinlat coslat ; 1 0 0 ; 0 coslat sinlat ] ;

% Transformation matrix from cylindrical to Cartesian basis
Rcarcyl = [ [coslon -sinlon 0] ; [sinlon coslon 0] ; [0 0 1] ] ;

u = (Rcarcyl * Rcylloc) * uloc;
