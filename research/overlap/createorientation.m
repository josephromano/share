function orient = createorientation(azDeg, altDeg);
%  CREATEORIENTATION -- Create orientation structure
%                       (converting angles from degrees to radians)
%
%  createorientation(azDeg, altDeg) creates a structure with fields
%        az: azimuth in radians East (clockwise) of North
%       alt: altitude (tilt) angle in radians above the local tangent plane
%  From the azimuth and altitude (in degrees).  If the second argument
%  is omitted, the altitude angle is set to zero.
%
%  This function simply converts the input angles from degrees into
%  radians (setting altDeg to zero if it's not included in the
%  argument list) and then packs them into a structure.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%  $Id: createorientation.m,v 1.3 2006/04/14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 1)
  altDeg = 0;
end
orient = struct('az',azDeg*pi/180,'alt',altDeg*pi/180);
