function r = getcartesianposition(loc);
%  GETCARTESIANPOSITION -- Construct Cartesian position vector from
%                          geographic coordinates
%
%  getcartesianposition(loc) Calculates the Cartesian position vector r 
%  (in meters) corresponding to a point on the Earth defined by the
%  location structure loc.  
%
%  The output is the position vector in units of meters in a Cartesian
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
%  The position is calculated assuming the latitude, longitude and
%  height are ellipsoidal coordinates on an ellipsoid of semimajor axis
%  6.378137e+6 meters and semiminor axis 6.356752314e+6.  These
%  coordinates are described in much greater detail in 
%  http://www.ligo.caltech.edu/docs/T/T980044-10.pdf
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%  $Id: getcartesianposition.m,v 1.3 2006/04/14 21:05:42 whelan Exp $
% 
%  See also CREATELOCATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is based on the WGS-84 ellipsoidal model of the Earth
a2 = 6.378137e+6 ^ 2; % square of semimajor axis (equatorial radius)
b2 = 6.356752314e+6 ^ 2; % square of semiminor axis (polar radius)
coslat = cos(loc.lat);
sinlat = sin(loc.lat);
R = sqrt(a2*coslat^2 + b2*sinlat^2);
R1 = (a2/R + loc.height)*coslat;
R2 = (b2/R + loc.height)*sinlat;

r = transpose([R1*cos(loc.lon) R1*sin(loc.lon) R2]);
