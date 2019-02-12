classdef Parameters
   % The definition of general information in the analysis 
   
   % Physical constants
   properties (Constant)
      kB = 1.38064852e-23
      me = 9.10938356e-31 % electron mass, [kg]
      mp = 1.6726219e-27  % proton mass, [kg]
      mi = 14             % average ion mass [amu]
      q  = 1.602e-19      % charge [C]
      Rg = 2634*1e3       % radius of Ganymede, [m]
   end
   
   % 
   properties (Constant)
      % Directory
      Dir        char = '~/Ganymede/MOP2018/runG8_PIC_1200s/EnergeticFlux/run_Analysis_largeDomain'
      % Follow B or anti-B direction
      %=========================================================
      % NOTE: user must set this and fnames, Region accordingly!
      %=========================================================
      Hemisphere char{mustBeMember(Hemisphere,{'north','south'})} = 'north'
      % Shell output containing surface B field info
      fnameSurf  char = 'shl_var_1_t00000557_n00250489.out'
      % Field info from PIC
      fnameField char = '3d_fluid_region0_0_t00000557_n00010710.out'
      % Electron
      fnameE     char = 'cut_particles0_region0_1_t00000557_n00010710.out'
      % Ion
      fnameI     char = 'cut_particles1_region0_2_t00000557_n00010710.out'
      %fnameI     char = 'cut_particles1_region0_4_t00000557_n00010710.out';
      % Topology info from MHD
      fnameGM    char = 'box_var_2_t00000557_n00250489.out'
      % Particle Species
      Species    char{mustBeMember(Species,{'electron','ion'})}= 'ion' 
      
      % Region of interest
      Region     = [-1.2 -1.125 -1.5 1.5 0.5 1.6];
      %Region     = [-1.2 -1.125 -1.5 1.5 -0.5 -1.6];
      
      % Preallocation size for particles
      ncountmax  double {mustBeInteger} = 1071080
      % Number of bins
      bins       double {mustBeInteger,mustBePositive} = 30
      
      % Index
      x_ = 1
      y_ = 2
      z_ = 3
      vx_= 4
      vy_= 5
      vz_= 6
      w_ = 7
   end
   
   % Unit conversion factors from runlog
   properties (Constant)
      Si2NoRho = 917591263.419815
      Si2NoV   = 2.5e-07
      Si2NoB   = 23.947746024154
      Si2NoE   = 5.9869365060385e-06
      Si2NoP   = 5.73494539637384e-05
      Si2NoJ   = 2.3947746024154e-06
      Si2NoL   = 1
      
      % Normalized to SI mass conversion
      No2SiMass = 1/(Parameters.Si2NoRho*Parameters.Si2NoL^3)
   end
   
   
end
