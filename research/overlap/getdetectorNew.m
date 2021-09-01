function [r, u, v, T] = getdetectorNew(site)
%
% GETDETECTORNEW -- a modified version of getdetector that returns
% detector geometry structure for lisa spacecraft vertices
%
%  The output is in the form of a structure with the fields
%      r: [3x1 double] %  position vector (in units of meters)
%                         in Cartesian coordinates
%      u: [3x1 double] %  unit vector along x-arm of detector
%                         in Cartesian coordinates
%      v: [3x1 double] %  unit vector along y-arm of detector
%                         in Cartesian coordinates
%      T:              %  arm length measured in light propagation
%                         time (in units of seconds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458;   % speed of light (m/s)

switch site

  case 'default'
    r = [0;0;0];
    u = [1;0;0];
    v = [0;1;0];
    T = 1;

  case 'lisaX'
    r = [0;0;0];
    u = [1;0;0];
    v = [cos(pi/3);sin(pi/3);0];
    T = 5e9/c;

  case 'lisaY'
    r = [5e9;0;0];
    u = [-cos(pi/3);sin(pi/3);0];
    v = [-1;0;0];
    T = 5e9/c;

  case 'lisaZ'
    r = [(5e9)*cos(pi/3);(5e9)*sin(pi/3);0];
    u = [-cos(pi/3);-sin(pi/3);0];
    v = [ cos(pi/3);-sin(pi/3);0];
    T = 5e9/c;

  case 'BBO1'
    r = [0;0;0];
    T = 5e9/c;
    u = [cos(pi/3);sin(pi/3);0];
    v = [-cos(pi/3);sin(pi/3);0];

  case 'BBO2'
    r = [0;(5e9)*2/sqrt(3);0];
    T = 5e9/c;
    u = [-cos(pi/3);-sin(pi/3);0];
    v = [cos(pi/3);-sin(pi/3);0];

  case 'L1'
    % LIGO Livinston 4km (Livingston, Louisiana, USA)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation( 30+(33+46.4196/60)/60 , ...
			  - (90+(46+27.2654/60)/60) , ...
   			  -6.574 );
    xarm = createorientation(180+72.2835, -3.121e-4*180/pi);
    yarm = createorientation(180-17.7165, -6.107e-4*180/pi);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = 3995.08/c;

  case 'H1'
    % LIGO Hanford 4km (Hanford, Washington, USA)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation( 46+(27+18.528/60)/60 , ...
			  - (119+(24+27.5657/60)/60) , ...
			  142.554);
    xarm = createorientation(-35.9994, -6.195e-4*180/pi);
    yarm = createorientation(180+54.0006, 1.25e-5*180/pi);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = 3995.08/c;

  case 'H2'
    % LIGO Hanford 2km (Hanford, Washington, USA)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation( 46+(27+18.528/60)/60 , ...
			  - (119+(24+27.5657/60)/60) , ...
			  142.554);
    xarm = createorientation(-35.9994, -6.195e-4*180/pi);
    yarm = createorientation(180+54.0006, 1.25e-5*180/pi);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = (3995.08/2)/c;

  case 'VIRGO'
    % VIRGO (Cascina/Pisa, Italy)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation(43 + (37 + 53.0921/60)/60 , ...
			  10 + (30 + 16.1878/60)/60 , ...
			  51.884);
    xarm = createorientation(90-70.5674);
    yarm = createorientation(90-160.5674);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = 3000/c; % EXACT NUMBER MIGHT BE DIFFERENT!!!

  case 'GEO600'
    % GEO-600 (Hannover, Germany)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation(52 + (14 + 42.528/60)/60 , ...
			  9 + (48 + 25.894/60)/60 , ...
			  114.425);
    xarm = createorientation(90-21.6117);
    yarm = createorientation(90-115.9431);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = 600/c; % EXACT NUMBER MIGHT BE DIFFERENT!!!

  case 'K1'
    % LCGT, Japan from Schutz, CQG 28 125023 (2011)

    loc  = createlocation(36 + (15 + 0/60)/60 , ...
                          137 + (10 + 48/60)/60 , ...
                          90);
    xarm = createorientation(65);
    yarm = createorientation(-25);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = 4000/c; % EXACT NUMBER MIGHT BE DIFFERENT!!!

  case 'I1'
    % LIGO-India from Schutz, CQG 28 125023 (2011)

    loc  = createlocation(19 + (5 + 47/60)/60 , ...
                          74 + (2 + 59/60)/60 , ...
                          90);
    xarm = createorientation(315);
    yarm = createorientation(225);
    r = getcartesianposition(loc);
    u = getcartesiandirection(xarm, loc);
    v = getcartesiandirection(yarm, loc);
    T = 4000/c; % EXACT NUMBER MIGHT BE DIFFERENT!!!

    % from sanjit, below
    %Det.lambda = [77 00 00]*[1; 1/60; 1/3600];
    %Det.phi = [13 00 00]*[1; 1/60; 1/3600];
    %Det.psi = [24 12 00]*[1; 1/60; 1/3600];

    %theta = (90-Det.phi)*pi/180;
    %phi   = (Det.lambda)*pi/180;
    %psi   = (Det.psi*pi)/180;

    %thetaCap = ...
    %    [ -cos(theta)*cos(phi); -cos(theta)*sin(phi); sin(theta) ];
    %phiCap = [ -sin(phi); cos(phi); 0.0 ];

    %Det.X    =  cos(psi)*thetaCap + sin(psi)*phiCap;
    %Det.Y    = -sin(psi)*thetaCap + cos(psi)*phiCap;
    %Det.Z    =  cross(Det.X,Det.Y);

    %----- Response matrix.
    %Det.d = 0.5*(Det.X*Det.X' - Det.Y*Det.Y') ;
    %Det.name = 'INDIGO';
    %Det.n = [];

  otherwise
    error('unrecognized detector site');

end

return
