classdef Parameters
   % The definition of general information in the analysis 
   
   % Physical constants
   properties (Constant)
      kb = 1.38064852e-23;
      m  = 1; % mass
      me = 9.10938356e-31; % electron mass, [kg]
      mp = 1.6726219e-27;  % proton mass, [kg]
      mi  = 14; % average ion mass [amu]
      q  = 1; % charge
      Rg = 2634*1e3; % radius of Ganymede, [m]
   end
   
   % 
   properties (Constant)
      % Directory
      Dir        char = '~/Documents/research/Ganymede/data/EnergeticFlux'
      % Shell output containing surface B field info
      fnameSurf  char = 'shl_var_1_t00000557_n00250489.out'
      % Field info from PIC
      fnameField char = '3d_fluid_region0_0_t00000557_n00010710.out'
      % Electron
      fnameE     char = 'cut_particles0_region0_1_t00000557_n00010710.out'
      % Ion
      fnameI     char = 'cut_particles1_region0_2_t00000557_n00010710.out'
      % Topology info from MHD
      fnameGM    char = 'box_var_2_t00000557_n00250489.out'
      % Particle Species
      Species    char{mustBeMember(Species,{'electron','ion'})}= 'ion' 
      
      % Region of interest
      Region     = [-1.2 -1.125 -2 2 0.5 2];
      % Preallocation size for particles
      ncountmax  double {mustBeInteger} = 1071080
   end

   
end