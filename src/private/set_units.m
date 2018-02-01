function set_units( filehead,type,varargin )
%set_units Set the units for the output files
%   If type is given as                                                            
%   'SI', 'CGS', 'NORMALIZED', 'PIC', 'PLANETARY', or 'SOLAR',                     
%   set typeunit = type otherwise try to guess from the fileheader.                
%
%   Based on typeunit set units for distance (xSI), time (tSI),                    
%   density (rhoSI), pressure (pSI), magnetic field (bSI)                          
%   and current density (jSI) in SI units.                                         
%   Distance unit (rplanet or rstar), ion and electron mass in amu                 
%   can be set with optional distunit, Mion and Melectron.                         
%
%   Also calculate convenient constants ti0, cs0 ... for typical formulas.         
%   See file "defaults" for definitions and usage 
%--------------------------------------------------------------------------

ndim      = filehead.ndim;
headline  = filehead.headline;
neqpar    = filehead.neqpar;
nw        = filehead.nw;
eqpar     = filehead.eqpar;
variables = filehead.variables;

% I can use a classdef to define all these constants if necessary.
% See my +constant directory for reference.
mu0SI = 4*pi*1e-7;      % H/m
cSI   = 2.9978e8;       % speed of light
mpSI  = 1.6726e-27;     % kg
eSI   = 1.602e-19;      % elementary charge, C
AuSI  = 149597870700;   % m
RsSI  = 6.957e8;        % m
Mi    = 1.0;            % Ion mass in amu
Me    = 1.0/1836.15;    % Electron mass in amu
gamma = 5/3;            % Adiabatic index for first fluid
gammae= 5/3;            % Adiabatic index for electrons
kbSI  = physconst('Boltzmann');
e0SI = 8.8542e-12;  % F/m

p = inputParser;
expectedtypes = {'SI','CGS','NORMALIZED','PIC','PLANETARY','SOLAR',' '};
defaultdistunit = [];      % default distance unit
defaultMion = [];          % default ion mass
defaultMelectron = [];     % default electron mass

addRequired(p,'filehead',@isstruct)
addRequired(p,'type',@(x) any(validatestring(x,expectedtypes)))

addOptional(p,'distunit',defaultdistunit,@isnumeric);
addParameter(p,'Mion',defaultMion,@isnumeric);
addParameter(p,'Melectron',defaultMelectron,@isnumeric);

parse(p,filehead,type,varargin{:});

% This part is used to guess the units.
% To be honest, I don`t understand the logic here. For example,
% nPa and m/s may appear in the same headline?
if ~strcmp(p.Results.type,' ')
   typeunit = upper(p.Results.type);
elseif contains(filehead.headline,'PIC')
   typeunit = 'PIC';
elseif contains(filehead.headline,' AU ')
   typeunit = 'OUTERHELIO';
elseif contains(filehead.headline,'kg/m3') || ...
      contains(filehead.headline,'m/s')
   typeunit = 'SI';
elseif contains(filehead.headline,'nPa') || ...
      contains(filehead.headline,' nT ')
   typeunit = 'PLANETARY';
elseif contains(filehead.headline,'dyne') || ...
      contains(filehead.headline,' G')
   typeunit = 'SOLAR';
else
   typeunit = 'NORMALIZED';
end

switch string(typeunit)
   case 'SI'
      xSI   = 1.0;             % m
      tSI   = 1.0;             % s
      rhoSI = 1.0;             % kg/m^3
      uSI   = 1.0;             % m/s
      pSI   = 1.0;             % Pa
      bSI   = 1.0;             % T
      jSI   = 1.0;             % A/m^2
   case 'CGS'
      xSI   = 0.01;            % cm
      tSI   = 1.0;             % s
      rhoSI = 1000.0;          % g/cm^3
      uSI   = 0.01;            % cm/s
      pSI   = 0.1;             % dyne/cm^2
      bSI   = 1.0e-4;          % G
      jSI   = 10*cSI;          % Fr/s/cm^2
   case 'PIC'
      % Normalized PIC units                                                           
      xSI   = 1.0;             % cm
      tSI   = 1.0;             % s
      rhoSI = 1.0;             % g/cm^3
      uSI   = 1.0;             % cm/s
      pSI   = 1.0;             % dyne/cm^2
      bSI   = 1.0;             % G
      jSI   = 1.0;             % Fr/s/cm^2
      c0    = 1.0;             % speed of light always 1 for iPIC3D
   case 'NORMALIZED'
      xSI   = 1.0;             % distance unit in SI
      tSI   = 1.0;             % time unit in SI
      rhoSI = 1.0;             % density unit in SI
      uSI   = 1.0;             % velocity unit in SI
      pSI   = 1.0;             % pressure unit in SI
      bSI   = sqrt(mu0SI);     % magnetic unit in SI
      jSI   = 1/sqrt(mu0SI);   % current unit in SI
      c0    = 1.0;             % speed of light (for Boris correction)
   case 'PLANETARY'
      xSI   = physconst('EarthRadius'); % Earth radius (default planet)
      tSI   = 1.0;             % s
      rhoSI = mpSI*1e6;        % mp/cm^3
      uSI   = 1e3;             % km/s
      pSI   = 1e-9;            % nPa
      bSI   = 1e-9;            % nT
      jSI   = 1e-6;            % muA/m^2
      c0    = cSI/uSI;         % speed of light in velocity units
   case 'OUTERHELIO'
      xSI   = AuSI;            % AU
      tSI   = 1.0;             % s
      rhoSI = mpSI*1e6;        % mp/cm^3
      uSI   = 1e3;             % km/s
      pSI   = 1e-1;            % dyne/cm^2
      bSI   = 1e-9;            % nT
      jSI   = 1e-6;            % muA/m^2
      c0    = cSI/uSI;         % speed of light in velocity units
   case 'SOLAR'
      xSI   = RsSI;            % radius of the Sun
      tSI   = 1.0;             % s
      rhoSI = 1e3;             % g/cm^3
      uSI   = 1e3;             % km/s
      pSI   = 1e-1;            % dyne/cm^2
      bSI   = 1e-4;            % G
      jSI   = 1e-6;            % muA/m^2
      c0    = cSI/uSI;         % speed of light in velocity units
   otherwise
      error('invalid typeunit=%s',typeunit)
end

% Overwrite values if given by eqpar
for ieqpar = 1:neqpar
   switch string(variables{ndim+nw+ieqpar})
      case 'xSI' 
         xSI   = eqpar(ieqpar);
      case 'tSI'
         tSI   = eqpar(ieqpar);
      case 'uSI'
         uSI   = eqpar(ieqpar);
      case 'rhoSI' 
         rhoSI = eqpar(ieqpar);
      case 'mi'
         mi    = eqpar(ieqpar);
      case 'm1'
         m1    = eqpar(ieqpar); 
      case 'me'
         me    = eqpar(ieqpar);
      case 'qi'
         qi    = eqpar(ieqpar);
      case 'q1'
         q1    = eqpar(ieqpar);
      case 'qe'
         qe    = eqpar(ieqpar); 
      case 'g'
         gamma = eqpar(ieqpar);
      case 'g1'
         gamma = eqpar(ieqpar);  
      case 'ge'
         ge    = eqpar(ieqpar);
      case 'c'
         c     = eqpar(ieqpar); 
      case 'clight'
         clight= eqpar(ieqpar);
      case 'r'
         r     = eqpar(ieqpar);  
      case 'rbody'
         rbody = eqpar(ieqpar);  
   end
end

% Overwrite distance unit if given as an argument 
if ~isempty(p.Results.distunit)
   xSI = p.Results.distunit;
end

% Overwrite ion and electron masses if given as an argument
if ~isempty(p.Results.Mion);      Mi = p.Results.Mion; end;
if ~isempty(p.Results.Melectron); Me = p.Results.Melectron; end;

% Calculate convenient conversion factors                                              
if strcmp(typeunit,'NORMALIZED')
   ti0  = 1.0/Mi;            % T      = p/rho*Mi           = ti0*p/rho
   cs0  = 1.0;               % cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
   mu0A = 1.0;               % vA     = sqrt(b/rho)        = sqrt(bb/mu0A/rho)
   mu0  = 1.0;               % beta   = p/(bb/2)           = p/(bb/(2*mu0))
   uH0  = Mi;                % uH     = j/rho*Mi           = uH0*j/rho
   op0  = 1.0/Mi;            % omegap = sqrt(rho)/Mi       = op0*sqrt(rho)
   oc0  = 1.0/Mi;            % omegac = b/Mi               = oc0*b
   rg0  = sqrt(Mi);          % rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
   di0  = c0*Mi;             % di = c0/sqrt(rho)*Mi        = di0/sqrt(rho)
   ld0  = Mi;                % ld = sqrt(p)/(rho*c0)*Mi    = ld0*sqrt(p)/rho
elseif strcmp(typeunit,'PIC')
   ti0  = 1.0/Mi;            % T      = p/rho*Mi           = ti0*p/rho
   cs0  = 1.0;               % cs     = sqrt(gamma*p/rho)  = sqrt(gs*p/rho)
   mu0A = 4*pi;              % vA     = sqrt(b/(4*!pi*rho))= sqrt(bb/mu0A/rho)
   mu0  = 4*pi;              % beta   = p/(bb/(8*!pi))     = p/(bb/(2*mu0))
   uH0  = Mi;                % uH     = j/rho*Mi           = uH0*j/rho
   op0  = sqrt(4*pi)/Mi;     % omegap = sqrt(4*!pi*rho)/Mi = op0*sqrt(rho)
   oc0  = 1.0/Mi;            % omegac = b/Mi               = oc0*b
   rg0  = sqrt(Mi);          % rg = sqrt(p/rho)/b*sqrt(Mi) = rg0*sqrt(p/rho)/b
   di0  = 1.0/sqrt(4*pi);    % di = 1/sqrt(4*!pi*rho)*Mi   = di0/sqrt(rho)
   ld0  = 1.0/sqrt(4*pi);    % ld = sqrt(p/(4*!pi))/rho*Mi = ld0*sqrt(p)/rho
else
   qom  = eSI/(Mi*mpSI); moq = 1/qom;
   ti0  = mpSI/kbSI*pSI/rhoSI*Mi;       % T[K]=p/(nk) = ti0*p/rho
   cs0  = pSI/rhoSI/uSI^2;              % cs          = sqrt(gs*p/rho)
   mu0A = uSI^2*mu0SI*rhoSI*bSI^(-2);   % vA          = sqrt(bb/(mu0A*rho))
   mu0  = mu0SI*pSI*bSI^(-2);           % beta        = p/(bb/(2*mu0))
   uH0  = moq*jSI/rhoSI/uSI;            % uH=j/(ne)   = uH0*j/rho
   op0  = qom*sqrt(rhoSI/e0SI)*tSI;     % omegap      = op0*sqrt(rho)
   oc0  = qom*bSI*tSI;                  % omegac      = oc0*b
   rg0  = moq*sqrt(pSI/rhoSI)/bSI/xSI/sqrt(Mi); % rg     = rg0*sqrt(p/rho)/b
   di0  = cSI/(op0/tSI)/xSI;                    % di=c/omegap = di0/sqrt(rho)
   ld0  = moq*sqrt(pSI)/rhoSI/xSI;              % ld          = ld0*sqrt(p)/rho
end

end