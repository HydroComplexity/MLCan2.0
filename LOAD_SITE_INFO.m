%************************** Number of Species ****************************%
load './Temps/temp_variable.mat'...
    'num_species' 'Sim_species' 'Sim_species_con'
% num_species : Number of species
% Sim_species = 1: Multi species
% Sim_species_con: Which species you would like to run

if Sim_species == 1
    kspecies = num_species;
else % for the case you would like to simulate a selected species
    kspecies = 5;
end

%***************************** Set Years *********************************%
load './Temps/temp_variable.mat'...
    'working_forcings'
load(working_forcings,...
    'year_crop');
Each_year = unique(year_crop)';

%************************ Set Initial Condition **************************%
load './Temps/temp_variable.mat'...
    'root_init' 'root_init_litter'
Tsinit=root_init(:,3);                          % Initial soil temperature ['C]
Tslint=root_init_litter(1,3);                   % Initial litter temperature ['C]
volliqinit=root_init(:,2);                      % Initial soil moisture [-]
volliqliinit=root_init_litter(1,2);             % Initial litter moisture [-]

%************ Site Independent Constant and Fixed Values *****************%
% CONVERSION FACTORS
CONSTANTS.umoltoWm2 = (2.17*10^5) / 10^6;       % Radiation conversion
CONSTANTS.Wm2toumol = 1/CONSTANTS.umoltoWm2;    % Radiation conversion
CONSTANTS.mmH2OtoMPa = 9.8066e-06;              % Pressure conversion

% PHYSICAL CONSTANTS
CONSTANTS.R = 8.314;                            % [J mol^-1 K^-1]
CONSTANTS.R_kJ = CONSTANTS.R/1000;              % [kJ mol^-1 K^-1]
CONSTANTS.Lv = 44000;                           % latent heat of vaporization [J/mol]
CONSTANTS.Lv_kg = 2260*1000;                    % Lv in [J / kg]
CONSTANTS.Lf = 6000;                            % latent heat of vaporization [J/mol]
CONSTANTS.Lf_kg = 3.337e5;                      % Lf in [J / kg]

CONSTANTS.cp_mol = 29.3;                        % specific heat of air at constant pressure [J/mol/K]
CONSTANTS.cp_JkgK = 1012;                       % [J/kg/K]
CONSTANTS.boltz = 5.6697 * 10^-8;               % Stefan-Boltzmann constant [W/m^2/K^4]
CONSTANTS.vonk = 0.41;                          % von Karman constant
% The K?m? constant is often used in turbulence modeling,
% for instance in boundary-layer meteorology to calculate fluxes
% of momentum, heat and moisture from the atmosphere to the land surface.
% It is considered to be a universal (? = 0.41).

CONSTANTS.rho_dry_air = 1.2923;                 % [kg / m^3]
CONSTANTS.grav = 9.8;                           % [m / s^2]

CONSTANTS.Vw = 18;                              % Mols to kg conversion

CONSTANTS.cpl_JkgK = 1800;                      % specific heat capacity of leaves[J/kg/K]
CONSTANTS.rho_leaf = 600;                       % Density of Leaves [kg / m^3]

CONSTANTS.fPAR = 0.5;                           % [-] fraction incoming shortwave as PAR

load './Temps/temp_variable.mat'...
    'days_step' 'mins_step' 'hours_step'
CONSTANTS.timestep = days_step*24*60+hours_step*60+mins_step; % [minutes]
CONSTANTS.dtime = CONSTANTS.timestep*60;        % [s]

% ROOT CUT INFORMATION
% wilting point
CONSTANTS.wilpoint = -1.5;                      % Wilting Point 1.5 MPa
% embolism
CONSTANTS.PLCa = 2.146;                         % Parameter a embolism equation
CONSTANTS.PLCb = -1.238;                        % Parameter b embolism equation
% Layer to cut
load './Temps/temp_variable.mat'...
    'set_para_root1'
CONSTANTS.nlc = set_para_root1{4,2};            % Number of layer to cut
VARIABLES.comroot = 0;                          % Initialize to zero the cut of roots

%*************************************************************************%
%********************** Site Specific Parameters *************************%
%*************************************************************************%
load './Temps/temp_variable.mat'...
    'num_LAD1'
if isstr(num_LAD1) == 1
    PARAMS.CanStruc.nl_can=str2num(num_LAD1);   % # canopy layers
else
    PARAMS.CanStruc.nl_can=(num_LAD1);          % # canopy layers
end

%*************************************************************************%
%                               Canopy                                    %
%*************************************************************************%
%********************* Independent of vegetation *************************%
% CANOPY STRUCTURE
load './Temps/temp_variable.mat'...
    'para_canopy_crop_fixed'
if isempty(cell2mat(para_canopy_crop_fixed(3,2)))
    PARAMS.CanStruc.hcan = cell2mat(para_canopy_crop_fixed(3,3));     % canopy height [m]
else
    PARAMS.CanStruc.hcan = cell2mat(para_canopy_crop_fixed(3,2));
end
if isempty(cell2mat(para_canopy_crop_fixed(4,2)))
    PARAMS.CanStruc.hobs = cell2mat(para_canopy_crop_fixed(4,3));     % observation height [m]
else
    PARAMS.CanStruc.hobs = cell2mat(para_canopy_crop_fixed(4,2));
end
if isempty(cell2mat(para_canopy_crop_fixed(5,2)))
    PARAMS.CanStruc.z0 = cell2mat(para_canopy_crop_fixed(5,3));       % canopy roughness length [m]
else
    PARAMS.CanStruc.z0 = cell2mat(para_canopy_crop_fixed(5,2));
end
PARAMS.CanStruc.d0 = 2/3 * PARAMS.CanStruc.hcan;                      % canopy displacement height [m]
PARAMS.CanStruc.VAratio = 3/1000;                                     % Ratio of volume to area [m]

% RADIATION
load './Temps/temp_variable.mat'...
    'para_radiation'
if isempty(cell2mat(para_radiation(1,2)))
    PARAMS.Rad.transmiss = cell2mat(para_radiation(1,3));             % atmospheric transmissivity
else
    PARAMS.Rad.transmiss = cell2mat(para_radiation(1,2));
end
if isempty(cell2mat(para_radiation(2,2)))
    PARAMS.Rad.epsv = cell2mat(para_radiation(2,3));                  % vegetation emissivity
else
    PARAMS.Rad.epsv = cell2mat(para_radiation(2,2));
end
if isempty(cell2mat(para_radiation(3,2)))
    PARAMS.Rad.epss = cell2mat(para_radiation(3,3));                  % soil emissivity
else
    PARAMS.Rad.epss = cell2mat(para_radiation(3,2));
end
if isempty(cell2mat(para_radiation(4,2)))
    PARAMS.Rad.epsa = cell2mat(para_radiation(4,3));                  % atmospheric emissivity
else
    PARAMS.Rad.epsa = cell2mat(para_radiation(4,2));
end
if isempty(cell2mat(para_radiation(5,2)))
    PARAMS.Rad.xx = cell2mat(para_radiation(5,3));                    % leaf angle dist param (Spherical is a good assumption (Shade and Golstein 2002)
else
    PARAMS.Rad.xx = cell2mat(para_radiation(5,2));
end
if isempty(cell2mat(para_radiation(6,2)))
    PARAMS.Rad.clump = cell2mat(para_radiation(6,3));                 % leaf clumping parameter
else
    PARAMS.Rad.clump = cell2mat(para_radiation(6,2));
end
if isempty(cell2mat(para_radiation(7,2)))
    PARAMS.Rad.Kdf = cell2mat(para_radiation(7,3));                   % extinction coeff for diffuse around 0.7 according to figure 15.2 (Campbell and Norman 1998)
else
    PARAMS.Rad.Kdf = cell2mat(para_radiation(7,2));
end

if isempty(cell2mat(para_radiation(8,2)))
    PARAMS.Rad.absorp_PAR = cell2mat(para_radiation(8,3));            % leaf absorptivity to PAR
else
    PARAMS.Rad.absorp_PAR = cell2mat(para_radiation(8,2));
end
if isempty(cell2mat(para_radiation(9,2)))
    PARAMS.Rad.absorp_NIR = cell2mat(para_radiation(9,3));            % leaf absorptivity to NIR
else
    PARAMS.Rad.absorp_NIR = cell2mat(para_radiation(9,2));
end
if isempty(cell2mat(para_radiation(10,2)))
    PARAMS.Rad.refl_PAR = cell2mat(para_radiation(10,3));             % PAR reflection coeff
else
    PARAMS.Rad.refl_PAR = cell2mat(para_radiation(10,2));
end
if isempty(cell2mat(para_radiation(11,2)))
    PARAMS.Rad.refl_NIR = cell2mat(para_radiation(11,3));             % NIR reflection coeff
else
    PARAMS.Rad.refl_NIR = cell2mat(para_radiation(11,2));
end
if isempty(cell2mat(para_radiation(12,2)))
    PARAMS.Rad.refl_soil = cell2mat(para_radiation(12,3));            % soil reflection coeff
else
    PARAMS.Rad.refl_soil = cell2mat(para_radiation(12,2));
end
PARAMS.Rad.trans_PAR = 1 - PARAMS.Rad.absorp_PAR - PARAMS.Rad.refl_PAR;
PARAMS.Rad.trans_NIR = 1 - PARAMS.Rad.absorp_NIR - PARAMS.Rad.refl_NIR;

if isempty(cell2mat(para_radiation(13,2)))
    PARAMS.Rad.refl_snow = cell2mat(para_radiation(13,3));            % Reflectivity of new snow.
else
    PARAMS.Rad.refl_snow = cell2mat(para_radiation(13,2));
end
if isempty(cell2mat(para_radiation(14,2)))
    PARAMS.Rad.refl_snow_old = cell2mat(para_radiation(14,3));        % How much reflectivity of snow decreases with age (per 12 days)
else
    PARAMS.Rad.refl_snow_old = cell2mat(para_radiation(14,2));
end

PARAMS.LWmethod = 0;                                                  % Method to compute in LW_ATTENUATION FUNCITON LW 0. Initial Darren. 1. Corrected
PARAMS.retainLW = 0;                                                  % Factor to compute LW, what proportion of LW radiation emmited by a canopy

% layer, remains there instead of leaving. Only if LWmethod =1
PARAMS.LWcom = 3;                                                     % Computation of longwave incoming 1.DATA 2.Boltzman 3.Boltzman Corrected

% RESPIRATION
load './Temps/temp_variable.mat'...
    'para_respiration' 'para_microenvironment'
if isempty(cell2mat(para_respiration(1,2)))
    PARAMS.Resp.Ro = cell2mat(para_respiration(1,3));                 % [umol / m^2 / s]
else
    PARAMS.Resp.Ro = cell2mat(para_respiration(1,2));
end
if isempty(cell2mat(para_respiration(2,2)))
    PARAMS.Resp.Q10 = cell2mat(para_respiration(2,3));
else
    PARAMS.Resp.Q10 = cell2mat(para_respiration(2,2));
end

% Turbulence
if isempty(cell2mat(para_microenvironment(1,2)))
    PARAMS.MicroEnv.Cd = cell2mat(para_microenvironment(1,3));        % Drag coefficient [-]
else
    PARAMS.MicroEnv.Cd = cell2mat(para_microenvironment(1,2));
end

PARAMS.MicroEnv.alph = CONSTANTS.vonk/3;                              % mixing length parameter [-] / Katul et al (BLM, 2004, p. 84)


load './Temps/temp_variable.mat'...
    'para_photosynthesisC3_crop1'
if isempty(cell2mat(para_photosynthesisC3_crop1(5,2)))
    PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop1(5,3)); % Vertical Distribution of Photosynthetic Capacity
else
    PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop1(5,2));
end

%*************************************************************************%
%                               SOIL                                      %
%*************************************************************************%
load './Temps/temp_variable.mat'...
    'num_root1' 'set_root_func' 'para_soil'
if isstr(num_root1) == 1
    PARAMS.nl_soil=str2num(num_root1);
else
    PARAMS.nl_soil=num_root1;
end

PARAMS.Soil.scaleza = cell2mat(set_root_func(1,2));
PARAMS.Soil.scalezb = cell2mat(set_root_func(2,2));

% SOIL PROPERTIES
if isempty(cell2mat(para_soil(1,2)))
    VERTSTRUC.sand = cell2mat(para_soil(1,3));
else
    VERTSTRUC.sand = cell2mat(para_soil(1,2));
end
if isempty(cell2mat(para_soil(2,2)))
    VERTSTRUC.clay = cell2mat(para_soil(2,3));
else
    VERTSTRUC.clay = cell2mat(para_soil(2,2));
end
if isempty(cell2mat(para_soil(3,2)))
    PARAMS.Soil.Cd_soil = cell2mat(para_soil(3,3));                   % soil drag coefficient [-]
else
    PARAMS.Soil.Cd_soil = cell2mat(para_soil(3,2));
end
if isempty(cell2mat(para_soil(4,2)))
    PARAMS.Soil.z0 = cell2mat(para_soil(4,3));                        % soil surface roughness length [m]
else
    PARAMS.Soil.z0 = cell2mat(para_soil(4,2));
end

PARAMS.Soil.alphCN=0.5;                                               % Cranck Nicholson for Heat Solution


% CONSTANTS
PARAMS.Soil.smpmin = -1.e8;                                           % restriction for min of soil poten. [mm]
PARAMS.Soil.wimp   = 0.05;                                            % water impremeable if porosity less than wimp [-]

PARAMS.Soil.scalek  = 0.5;                                            % length scale for the exponential decrease in Ksat [m]


% HEAT CAPACITIES [J / kg / K)]
% Dry Air
PARAMS.Soil.HC_air = 1.00464 * 10^3;
% Water
PARAMS.Soil.HC_liq = 4.188 * 10^3;
% Ice
PARAMS.Soil.HC_ice = 2.11727 * 10^3;
% DENSITIES [kg / m^3]
PARAMS.Soil.rho_liq = 1000;
PARAMS.Soil.rho_ice = 917;
% THERMAL CONDUCTIVITIES  [W / m / K]
PARAMS.Soil.TK_liq = 0.6;
PARAMS.Soil.TK_ice = 2.29;
PARAMS.Soil.TK_air = 0.023;
% FREEZING TEMP OF FRESH WATER [K]
PARAMS.Soil.Tf = 273.16;

% FOR LITTER SNOW
PARAMS.Tcr_atm = 273;                                                 % Atmospheric temperature above which snow
PARAMS.Soil.HC_snow = 0.1;                                            % Snow holding capacity of liquid water
PARAMS.Soil.rhosn_min = 50;                                           % Snow compaction parameter [kg/m3]
PARAMS.Soil.fcomp = 0.6;                                              %
PARAMS.Soil.rhod = 150.00;                                            % Snow compaction parameter [kg/m3]

PARAMS.Soil.c1 = -1.4e4;
PARAMS.Soil.thetals = 0.9;                                            % value of soil moisture litter at saturation
PARAMS.Soil.thetafc = 0.025;                                          % Value of litter soil moisture at field capacity
PARAMS.Soil.thetatr = 0.18;                                           % Value of soil moisture after which
% rl becomes negligible
PARAMS.Soil.psill = 35.3;                                             % parameter to compute psi litter
PARAMS.Soil.bl = 2.42;                                                % parameter to compute psi litter
PARAMS.Soil.bdl = 42.5;                                               % Bulk density of litter [kg/m3]
PARAMS.Soil.rhowater = 1000;                                          % Density of liquid water [1000 kg/m3]
PARAMS.Soil.TK_litter = 0.15;                                         %  Litter Thermal Conductivity [W/m/k] or [J/s/m/k]
PARAMS.Soil.TD_litter = 5.7*10^(-7);                                  % Litter thermal diffusivity  [m/s]
PARAMS.Soil.slfactor = 0.1;                                           % Percentage above which is considered as snow for energy balance
PARAMS.Soil.VHC_litter = 0.3*10^(6);                                  % Volumetric Heat Capacity of Litter [J/m3/K]
PARAMS.Soil.km = 0.0001;%0.00004;                                     % parameter to compute drainage from litter
PARAMS.Soil.bm = 0.1;%2.3;                                            % parameter to compute drainage from litter
PARAMS.Soil.thetamin = 0.0001;                                        % value of soil moisture at which evaporation becomes negligible
PARAMS.Soil.ldif = 3.1*10^(-8);
PARAMS.Soil.sdif = 2.12*10^(-5);
PARAMS.Soil.kklitter = 1000;                                          % [1/cm] Radiation attenuation
PARAMS.Soil.kksnow = 0.8;                                             % [1/cm] Radiation attenuation

load './Temps/temp_variable.mat'...
    'litter_depth'
if SWITCHES.litter
    VARIABLES.SOIL.dzlit_m = litter_depth;                            % Litter thickness in [m]
else
    VARIABLES.SOIL.dzlit_m = 0.;                                      % Litter thickness in [m]
end

% ENTROPY COMPUTATION.
PARAMS.Entropy.c2 = 2.336;
PARAMS.Entropy.c3 = 0.26;

%*************************************************************************%
%********************* DEPENDENT OF VEGETATION TYPE **********************%
%*************************************************************************%
% kspecies = 1.  Species 1
% kspecies = 2.  Species 1+2
% kspecies = 3.  Species 1+2+3
% kspecies = 4.  Species 1+2+3+4
% kspecies = 5.  Defined species

if (kspecies == 1 || kspecies == 5)
    PARAMS.CanStruc.nspecies = 1;  % Number of species
    nspecies = 1;
elseif (kspecies == 2)
    PARAMS.CanStruc.nspecies = 2;  % Number of species
    nspecies = 2;
elseif (kspecies == 3)
    PARAMS.CanStruc.nspecies = 3;  % Number of species
    nspecies = 3;
elseif (kspecies == 4 )
    PARAMS.CanStruc.nspecies = 4;  % Number of species
    nspecies = 4;
end

%**************************** Species 1 **********************************%
ii=0;
if (kspecies == 1 || kspecies == 2 || kspecies == 3 || kspecies == 4 || (kspecies == 5 && Sim_species_con == 1 ))
    ii = 1;
elseif ((kspecies == 5 && Sim_species_con == 2) || (kspecies == 5 && Sim_species_con == 3) || (kspecies == 5 && Sim_species_con == 4))
    ii = 0;
end

if ii ~= 0
    load './Temps/temp_variable.mat'...
        'para_leaf_crop1'
    
    if isempty(cell2mat(para_leaf_crop1(1,2)))
        PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop1(1,3));  % multiplicative factor for Fc and LE
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop1(1,2))) == 1
            PARAMS.CanStruc.LEfact(ii) = str2num(cell2mat(para_leaf_crop1(1,2)));
        else
            PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop1(1,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop1(2,2)))
        PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop1(2,3));   % multiplicative factor for H and LW
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
        
    else
        if isstr(cell2mat(para_leaf_crop1(2,2))) == 1
            PARAMS.CanStruc.Hfact(ii) = str2num(cell2mat(para_leaf_crop1(2,2)));
        else
            PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop1(2,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop1(3,2)))
        PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop1(3,3));
    else
        if isstr(cell2mat(para_leaf_crop1(3,2))) == 1
            PARAMS.CanStruc.LWfact(ii) = str2num(cell2mat(para_leaf_crop1(3,2)));
        else
            PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop1(3,2));
        end
    end
    
    % CANOPY STRUCTURE
    if isempty(cell2mat(para_leaf_crop1(4,2)))
        PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop1(4,3)); % 1 = broad leaves, 2 = needles
    else
        if isstr(cell2mat(para_leaf_crop1(4,2))) == 1
            PARAMS.CanStruc.leaftype(ii) = str2num(cell2mat(para_leaf_crop1(4,2)));
        else
            PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop1(4,2));
        end
    end
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop1' 'para_leaf_crop11'
    if isempty(cell2mat(para_canopy_crop1(1,2)))
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop1(1,3));     % leaf width or needle diameter [m]
    else
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop1(1,2));
    end
    if isempty(cell2mat(para_canopy_crop1(2,2)))
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop1(2,3));     % shoot diameter for conifers
        % leaf width for broadleaved vegetation (= ld)
    else
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop1(2,2));
    end
    
    % PARAMS WEIBULL / Dongkook: You can ignore below - Just keep here
    PARAMS.CanStruc.beta(ii) = 6.6533;%2.1681;
    PARAMS.CanStruc.alpha(ii) = 0.8017;%0.5831;
    
    % Smax = maximum h2o storage capacity for foliage [mm/LAI unit]
    if isempty(cell2mat(para_leaf_crop11(1,2)))
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop11(1,3));
    else
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop11(1,2));
    end
    
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop_fixed'
    if isempty(cell2mat(para_canopy_crop_fixed(1,2)))
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,3));      % max fraction of canopy that can be wet
    else
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,2));
    end
    if isempty(cell2mat(para_canopy_crop_fixed(2,2)))
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,3)); % precip extinction coefficient
    else
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,2));
    end
    
    % PHOTOSYNTHESIS
    load './Temps/temp_variable.mat'...
        'fullpath_forcings' 'para_canopy_crop1' 'para_photosynthesisC3_crop1' 'para_photosynthesisC4_crop1' 'ph_type1'
    load (fullpath_forcings)
    
    PARAMS.Photosyn.ph_type(ii) = ph_type1;
    
    if PARAMS.Photosyn.ph_type(ii) == 1
        % C3 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC3_crop1(1,2)))
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop1(1,3));
        else
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop1(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop1(2,2)))
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop1(2,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop1(2,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop1(3,2)))
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop1(3,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop1(3,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop1(4,2)))
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop1(4,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop1(4,2))*ones(1,size(Ta_crop,1));
        end
        
        if isempty(cell2mat(para_photosynthesisC3_crop1(5,2)))
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop1(5,3)); % Vertical Distribution of Photosynthetic Capacity
        else
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop1(5,2));
        end
        
        % C4 Photosynthesis Parameters
        PARAMS.Photosyn.Vmax_C4(ii) = NaN;
        PARAMS.Photosyn.Rd_C4(ii) = NaN;
        PARAMS.Photosyn.Q10_C4(ii) = NaN;
        PARAMS.Photosyn.kk_C4(ii) = NaN;
        PARAMS.Photosyn.theta_C4(ii) = NaN;
        PARAMS.Photosyn.beta_C4(ii) = NaN;
        PARAMS.Photosyn.al_C4(ii) = NaN;
    else
        % C3 Photosynthesis Parameters
        PARAMS.Photosyn.beta_ph_C3(ii) = NaN;
        PARAMS.Photosyn.Vcmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Jmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Rd25{ii} = NaN;
        
        % C4 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC4_crop1(1,2)))
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop1(1,3));
        else
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop1(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop1(2,2)))
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop1(2,3));
        else
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop1(2,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop1(3,2)))
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop1(3,3));
        else
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop1(3,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop1(4,2)))
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop1(4,3));
        else
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop1(4,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop1(5,2)))
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop1(5,3));
        else
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop1(5,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop1(6,2)))
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop1(6,3));
        else
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop1(6,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop1(7,2)))
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop1(7,3));
        else
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop1(7,2));
        end
        
        if isempty(cell2mat(para_photosynthesisC4_crop1(8,2)))
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC4_crop1(8,3)); % Vertical Distribution of Photosynthetic Capacity
        else
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC4_crop1(8,2));
        end
    end
    
    
    PARAMS.Photosyn.Oi(ii) = 210;                                           % Intercellular oxygen concentration [mmol / mol]
    
    
    % STOMATAL CONDUCTANCE
    load './Temps/temp_variable.mat'...
        'para_conductance_crop1'
    
    % Ball-Berry
    if isempty(cell2mat(para_conductance_crop1(1,2)))
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop1(1,3)); % slope parameter in BB model [-] (Leakey: m=10.6 (ambient); m=10.9 (elevated))
    else
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop1(1,2));
    end
    if isempty(cell2mat(para_conductance_crop1(2,2)))
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop1(2,3));   % intercept parameter in BB model [mol/m^2/s] (Leakey: b=0.008 (ambient); b=0.007 (elevated))
    else
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop1(2,2));
    end
    
    
    % (Tuzet et al, PCE 2003)
    if isempty(cell2mat(para_conductance_crop1(3,2)))
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop1(3,3));    % sensitivity parameter for initial decrease in leaf potential [-]
    else
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop1(3,2));
    end
    if isempty(cell2mat(para_conductance_crop1(4,2)))
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop1(4,3));  % leaf potential at which half of the hydraulic conductance is lost [MPa]
    else
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop1(4,2));
    end
    
    
    % Conductivities
    if isempty(cell2mat(para_conductance_crop1(5,2)))
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop1(5,3));      % radial conductivity of the root system [s^-1]
    else
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop1(5,2));
    end
    if isempty(cell2mat(para_conductance_crop1(6,2)))
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop1(6,3));      % axial specific conductivity of the root system [mm/s]
    else
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop1(6,2));
    end
    if isempty(cell2mat(para_conductance_crop1(7,2)))
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop1(7,3));
    else
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop1(7,2));
    end
    
    if isempty(cell2mat(para_conductance_crop1(8,2)))
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop1(8,3));
    else
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop1(8,2));
    end
    
    % ROOT SOIL
    % HR
    load './Temps/temp_variable.mat'...
        'HR'
    
    SWITCHES.HR_on(ii) = HR;                                          % Hydraulic Redistribution 1. Yes, 0. No
    
    % root structure
    load './Temps/temp_variable.mat'...
        'set_para_root1'
    PARAMS.Soil.z50f(ii) = cell2mat(set_para_root1(2,2));             % From Ameriflux site description
    PARAMS.Soil.z95f(ii) = cell2mat(set_para_root1(3,2));             % From Ameriflux site description
    PARAMS.Soil.z50t(ii) = cell2mat(set_para_root1(2,2));
    PARAMS.Soil.z95t(ii) = cell2mat(set_para_root1(3,2));
    PARAMS.Soil.maxrootdepth(ii) = cell2mat(set_para_root1(1,2));
    
end


%**************************** Species 2 **********************************%
ii=0;
if (kspecies == 5 && Sim_species_con == 2)
    ii = 1;
elseif (kspecies == 2 || kspecies == 3 || kspecies == 4)
    ii = 2;
elseif ((kspecies == 5 && Sim_species_con == 1) || (kspecies == 5 && Sim_species_con == 3) || (kspecies == 5 && Sim_species_con == 4))
    ii = 0;
end

if ii~= 0
    load './Temps/temp_variable.mat'...
        'para_leaf_crop2'
    
    if isempty(cell2mat(para_leaf_crop2(1,2)))
        PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop2(1,3));  % multiplicative factor for Fc and LE
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop2(1,2))) == 1
            PARAMS.CanStruc.LEfact(ii) = str2num(cell2mat(para_leaf_crop2(1,2)));
        else
            PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop2(1,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop2(2,2)))
        PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop2(2,3));   % multiplicative factor for H and LW
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop2(2,2))) == 1
            PARAMS.CanStruc.Hfact(ii) = str2num(cell2mat(para_leaf_crop2(2,2)));
        else
            PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop2(2,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop2(3,2)))
        PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop2(3,3));
    else
        if isstr(cell2mat(para_leaf_crop2(3,2))) == 1
            PARAMS.CanStruc.LWfact(ii) = str2num(cell2mat(para_leaf_crop2(3,2)));
        else
            PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop2(3,2));
        end
    end
    
    % CANOPY STRUCTURE
    if isempty(cell2mat(para_leaf_crop2(4,2)))
        PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop2(4,3)); % 1 = broad leaves, 2 = needles
    else
        if isstr(cell2mat(para_leaf_crop2(4,2))) == 1
            PARAMS.CanStruc.leaftype(ii) = str2num(cell2mat(para_leaf_crop2(4,2)));
        else
            PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop2(4,2));
        end
    end
    
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop2' 'para_leaf_crop22'
    
    if isempty(cell2mat(para_canopy_crop2(1,2)))
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop2(1,3));    % leaf width or needle diameter [m]
    else
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop2(1,2));
    end
    if isempty(cell2mat(para_canopy_crop2(2,2)))
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop2(2,3));    % shoot diameter for conifers
        % leaf width for broadleaved vegetation (= ld)
    else
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop2(2,2));
    end
    
    
    % PARAMS WEIBULL / Dongkook: You can ignore below - Just keep here
    PARAMS.CanStruc.beta(ii) = 2.1681;%2.1681;
    PARAMS.CanStruc.alpha(ii) = 0.5831;
    
    % Smax = maximum h2o storage capacity for foliage [mm/LAI unit]
    if isempty(cell2mat(para_leaf_crop22(1,2)))
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop22(1,3));
    else
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop22(1,2));
    end
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop_fixed'
    
    if isempty(cell2mat(para_canopy_crop_fixed(1,2)))
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,3)); % max fraction of canopy that can be wet
    else
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,2));
    end
    if isempty(cell2mat(para_canopy_crop_fixed(2,2)))
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,3)); % precip extinction coefficient
    else
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,2));
    end
    
    % PHOTOSYNTHESIS
    load './Temps/temp_variable.mat'...
        'fullpath_forcings' 'para_canopy_crop2' 'para_photosynthesisC3_crop2' 'para_photosynthesisC4_crop2' 'ph_type2'
    load (fullpath_forcings)
    
    PARAMS.Photosyn.ph_type(ii) = ph_type2;
    
    if PARAMS.Photosyn.ph_type(ii) == 1
        % C3 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC3_crop2(1,2)))
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop2(1,3));
        else
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop2(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop2(2,2)))
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop2(2,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop2(2,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop2(3,2)))
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop2(3,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop2(3,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop2(4,2)))
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop2(4,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop2(4,2))*ones(1,size(Ta_crop,1));
        end
        
        % C4 Photosynthesis Parameters
        PARAMS.Photosyn.Vmax_C4(ii) = NaN;
        PARAMS.Photosyn.Rd_C4(ii) = NaN;
        PARAMS.Photosyn.Q10_C4(ii) = NaN;
        PARAMS.Photosyn.kk_C4(ii) = NaN;
        PARAMS.Photosyn.theta_C4(ii) = NaN;
        PARAMS.Photosyn.beta_C4(ii) = NaN;
        PARAMS.Photosyn.al_C4(ii) = NaN;
    else
        % C3 Photosynthesis Parameters
        PARAMS.Photosyn.beta_ph_C3(ii) = NaN;
        PARAMS.Photosyn.Vcmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Jmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Rd25{ii} = NaN;
        
        % C4 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC4_crop2(1,2)))
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop2(1,3));
        else
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop2(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop2(2,2)))
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop2(2,3));
        else
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop2(2,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop2(3,2)))
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop2(3,3));
        else
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop2(3,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop2(4,2)))
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop2(4,3));
        else
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop2(4,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop2(5,2)))
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop2(5,3));
        else
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop2(5,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop2(6,2)))
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop2(6,3));
        else
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop2(6,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop2(7,2)))
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop2(7,3));
        else
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop2(7,2));
        end
    end
    
    PARAMS.Photosyn.Oi(ii) = 210;                                           % Intercellular oxygen concentration [mmol / mol]
    
    
    % STOMATAL CONDUCTANCE
    load './Temps/temp_variable.mat'...
        'para_conductance_crop2'
    
    % Ball-Berry
    if isempty(cell2mat(para_conductance_crop2(1,2)))
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop2(1,3)); % slope parameter in BB model [-] (Leakey: m=10.6 (ambient); m=10.9 (elevated))
    else
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop2(1,2));
    end
    if isempty(cell2mat(para_conductance_crop2(2,2)))
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop2(2,3));   % intercept parameter in BB model [mol/m^2/s] (Leakey: b=0.008 (ambient); b=0.007 (elevated))
    else
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop2(2,2));
    end
    
    % (Tuzet et al, PCE 2003)
    if isempty(cell2mat(para_conductance_crop2(3,2)))
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop2(3,3));    % sensitivity parameter for initial decrease in leaf potential [-]
    else
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop2(3,2));
    end
    if isempty(cell2mat(para_conductance_crop2(4,2)))
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop2(4,3));  % leaf potential at which half of the hydraulic conductance is lost [MPa]
    else
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop2(4,2));
    end
    
    % Conductivities
    if isempty(cell2mat(para_conductance_crop2(5,2)))
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop2(5,3));      % radial conductivity of the root system [s^-1]
    else
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop2(5,2));
    end
    if isempty(cell2mat(para_conductance_crop2(6,2)))
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop2(6,3));      % axial specific conductivity of the root system [mm/s]
    else
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop2(6,2));
    end
    if isempty(cell2mat(para_conductance_crop2(7,2)))
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop2(7,3));
    else
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop2(7,2));
    end
    
    if isempty(cell2mat(para_conductance_crop2(8,2)))
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop2(8,3));
    else
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop2(8,2));
    end
    
    % ROOT SOIL
    % HR
    load './Temps/temp_variable.mat'...
        'HR'
    SWITCHES.HR_on(ii) = HR;                                                % Hydraulic Redistribution 1. Yes, 0. No
    
    % root structure
    load './Temps/temp_variable.mat'...
        'set_para_root2'
    PARAMS.Soil.z50f(ii) = cell2mat(set_para_root2(2,2));                   % From Ameriflux site description
    PARAMS.Soil.z95f(ii) = cell2mat(set_para_root2(3,2));                   % From Ameriflux site description
    PARAMS.Soil.z50t(ii) = cell2mat(set_para_root2(2,2));
    PARAMS.Soil.z95t(ii) = cell2mat(set_para_root2(3,2));
    PARAMS.Soil.maxrootdepth(ii) = cell2mat(set_para_root2(1,2));
end


%**************************** Species 3 **********************************%
ii=0;
if (kspecies == 5 && Sim_species_con == 3)
    ii = 1;
elseif (kspecies == 3 || kspecies == 4)
    ii = 3;
elseif ((kspecies == 5 && Sim_species_con == 1) || (kspecies == 5 && Sim_species_con == 2) || (kspecies == 5 && Sim_species_con == 4))
    ii = 0;
end

if ii~=0
    load './Temps/temp_variable.mat'...
        'para_leaf_crop3'
    if isempty(cell2mat(para_leaf_crop3(1,2)))
        PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop3(1,3));  % multiplicative factor for Fc and LE
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop3(1,2))) == 1
            PARAMS.CanStruc.LEfact(ii) = str2num(cell2mat(para_leaf_crop3(1,2)));
        else
            PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop3(1,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop3(2,2)))
        PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop3(2,3));   % multiplicative factor for H and LW
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop3(2,2))) == 1
            PARAMS.CanStruc.Hfact(ii) = str2num(cell2mat(para_leaf_crop3(2,2)));
        else
            PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop3(2,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop3(3,2)))
        PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop3(3,3));
    else
        if isstr(cell2mat(para_leaf_crop3(3,2))) == 1
            PARAMS.CanStruc.LWfact(ii) = str2num(cell2mat(para_leaf_crop3(3,2)));
        else
            PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop3(3,2));
        end
    end
    
    % CANOPY STRUCTURE
    if isempty(cell2mat(para_leaf_crop3(4,2)))
        PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop3(4,3)); % 1 = broad leaves, 2 = needles
    else
        if isstr(cell2mat(para_leaf_crop3(4,2))) == 1
            PARAMS.CanStruc.leaftype(ii) = str2num(cell2mat(para_leaf_crop3(4,2)));
        else
            PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop3(4,2));
        end
    end
    
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop3' 'para_leaf_crop33'
    if isempty(cell2mat(para_canopy_crop3(1,2)))
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop3(1,3));    % leaf width or needle diameter [m]
    else
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop3(1,2));
    end
    if isempty(cell2mat(para_canopy_crop3(2,2)))
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop3(2,3));    % shoot diameter for conifers
        % leaf width for broadleaved vegetation (= ld)
    else
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop3(2,2));
    end
    
    
    % PARAMS WEIBULL / Dongkook: You can ignore below - Just keep here
    PARAMS.CanStruc.beta(ii) = 2.1681;
    PARAMS.CanStruc.alpha(ii) = 0.5831;
    
    % Smax = maximum h2o storage capacity for foliage [mm/LAI unit]
    if isempty(cell2mat(para_leaf_crop33(1,2)))
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop33(1,3));
    else
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop33(1,2));
    end
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop_fixed'
    
    if isempty(cell2mat(para_canopy_crop_fixed(1,2)))
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,3)); % max fraction of canopy that can be wet
    else
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,2));
    end
    if isempty(cell2mat(para_canopy_crop_fixed(2,2)))
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,3)); % precip extinction coefficient
    else
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,2));
    end
    
    % PHOTOSYNTHESIS
    load './Temps/temp_variable.mat'...
        'fullpath_forcings' 'para_canopy_crop3' 'para_photosynthesisC3_crop3' 'para_photosynthesisC4_crop3' 'ph_type3'
    load (fullpath_forcings)
    
    PARAMS.Photosyn.ph_type(ii) = ph_type3;
    
    if PARAMS.Photosyn.ph_type(ii) == 1
        % C3 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC3_crop3(1,2)))
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop3(1,3));
        else
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop3(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop3(2,2)))
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop3(2,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop3(2,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop3(3,2)))
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop3(3,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop3(3,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop3(4,2)))
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop3(4,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop3(4,2))*ones(1,size(Ta_crop,1));
        end
        
        % C4 Photosynthesis Parameters
        PARAMS.Photosyn.Vmax_C4(ii) = NaN;
        PARAMS.Photosyn.Rd_C4(ii) = NaN;
        PARAMS.Photosyn.Q10_C4(ii) = NaN;
        PARAMS.Photosyn.kk_C4(ii) = NaN;
        PARAMS.Photosyn.theta_C4(ii) = NaN;
        PARAMS.Photosyn.beta_C4(ii) = NaN;
        PARAMS.Photosyn.al_C4(ii) = NaN;
    else
        % C3 Photosynthesis Parameters
        PARAMS.Photosyn.beta_ph_C3(ii) = NaN;
        PARAMS.Photosyn.Vcmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Jmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Rd25{ii} = NaN;
        
        % C4 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC4_crop3(1,2)))
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop3(1,3));
        else
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop3(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop3(2,2)))
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop3(2,3));
        else
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop3(2,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop3(3,2)))
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop3(3,3));
        else
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop3(3,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop3(4,2)))
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop3(4,3));
        else
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop3(4,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop3(5,2)))
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop3(5,3));
        else
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop3(5,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop3(6,2)))
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop3(6,3));
        else
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop3(6,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop3(7,2)))
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop3(7,3));
        else
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop3(7,2));
        end
    end
    
    PARAMS.Photosyn.Oi(ii) = 210;                                     % Intercellular oxygen concentration [mmol / mol]
    
    
    % STOMATAL CONDUCTANCE
    load './Temps/temp_variable.mat'...
        'para_conductance_crop3'
    
    % Ball-Berry
    if isempty(cell2mat(para_conductance_crop3(1,2)))
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop3(1,3)); % slope parameter in BB model [-] (Leakey: m=10.6 (ambient); m=10.9 (elevated))
    else
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop3(1,2));
    end
    if isempty(cell2mat(para_conductance_crop3(2,2)))
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop3(2,3)); % intercept parameter in BB model [mol/m^2/s] (Leakey: b=0.008 (ambient); b=0.007 (elevated))
    else
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop3(2,2));
    end
    
    % (Tuzet et al, PCE 2003)
    if isempty(cell2mat(para_conductance_crop3(3,2)))
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop3(3,3)); % sensitivity parameter for initial decrease in leaf potential [-]
    else
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop3(3,2));
    end
    if isempty(cell2mat(para_conductance_crop3(4,2)))
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop3(4,3)); % leaf potential at which half of the hydraulic conductance is lost [MPa]
    else
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop3(4,2));
    end
    
    % Conductivities
    if isempty(cell2mat(para_conductance_crop3(5,2)))
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop3(5,3)); % radial conductivity of the root system [s^-1]
    else
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop3(5,2));
    end
    if isempty(cell2mat(para_conductance_crop3(6,2)))
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop3(6,3)); % axial specific conductivity of the root system [mm/s]
    else
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop3(6,2));
    end
    if isempty(cell2mat(para_conductance_crop3(7,2)))
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop3(7,3));
    else
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop3(7,2));
    end
    
    if isempty(cell2mat(para_conductance_crop3(8,2)))
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop3(8,3));
    else
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop3(8,2));
    end
    
    
    % ROOT SOIL
    % HR
    load './Temps/temp_variable.mat'...
        'HR'
    SWITCHES.HR_on(ii) = HR;                                          % Hydraulic Redistribution 1. Yes, 0. No
    
    % root structure
    load './Temps/temp_variable.mat'...
        'set_para_root3'
    PARAMS.Soil.z50f(ii) = cell2mat(set_para_root3(2,2));             % From Ameriflux site description
    PARAMS.Soil.z95f(ii) = cell2mat(set_para_root3(3,2));             % From Ameriflux site description
    PARAMS.Soil.z50t(ii) = cell2mat(set_para_root3(2,2));
    PARAMS.Soil.z95t(ii) = cell2mat(set_para_root3(3,2));
    PARAMS.Soil.maxrootdepth(ii) = cell2mat(set_para_root3(1,2));
end


%**************************** Species 4 **********************************%
ii=0;
if (kspecies == 5 && Sim_species_con == 4)
    ii = 1;
elseif (kspecies == 4)
    ii = 4;
elseif ((kspecies == 5 && Sim_species_con == 1) || (kspecies == 5 && Sim_species_con == 2) || (kspecies == 5 && Sim_species_con == 3))
    ii = 0;
end

if ii~=0
    
    
    load './Temps/temp_variable.mat'...
        'para_leaf_crop4'
    
    if isempty(cell2mat(para_leaf_crop4(1,2)))
        PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop4(1,3));  % multiplicative factor for Fc and LE
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop4(1,2))) == 1
            PARAMS.CanStruc.LEfact(ii) = str2num(cell2mat(para_leaf_crop4(1,2)));
        else
            PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop4(1,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop4(2,2)))
        PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop4(2,3));   % multiplicative factor for H and LW
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop4(2,2))) == 1
            PARAMS.CanStruc.Hfact(ii) = str2num(cell2mat(para_leaf_crop4(2,2)));
        else
            PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop4(2,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop4(3,2)))
        PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop4(3,3));
    else
        if isstr(cell2mat(para_leaf_crop4(3,2))) == 1
            PARAMS.CanStruc.LWfact(ii) = str2num(cell2mat(para_leaf_crop4(3,2)));
        else
            PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop4(3,2));
        end
    end
    
    % CANOPY STRUCTURE
    if isempty(cell2mat(para_leaf_crop4(4,2)))
        PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop4(4,3)); % 1 = broad leaves, 2 = needles
    else
        if isstr(cell2mat(para_leaf_crop4(4,2))) == 1
            PARAMS.CanStruc.leaftype(ii) = str2num(cell2mat(para_leaf_crop4(4,2)));
        else
            PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop4(4,2));
        end
    end
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop4' 'para_leaf_crop44'
    
    if isempty(cell2mat(para_canopy_crop4(1,2)))
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop4(1,3));    % leaf width or needle diameter [m
    else
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop4(1,2));
    end
    if isempty(cell2mat(para_canopy_crop4(2,2)))
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop4(2,3));    % shoot diameter for conifers
        % leaf width for broadleaved vegetation (= ld)
    else
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop4(2,2));
    end
    
    
    % PARAMS WEIBULL / Dongkook: You can ignore below - Just keep here
    PARAMS.CanStruc.beta(ii) = 2.1681;
    PARAMS.CanStruc.alpha(ii) = 0.5831;
    
    % Smax = maximum h2o storage capacity for foliage [mm/LAI unit]
    if isempty(cell2mat(para_leaf_crop44(1,2)))
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop44(1,3));
    else
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop44(1,2));
    end
    
    load './Temps/temp_variable.mat'...
        'para_canopy_crop_fixed'
    
    if isempty(cell2mat(para_canopy_crop_fixed(1,2)))
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,3)); % max fraction of canopy that can be wet
    else
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,2));
    end
    if isempty(cell2mat(para_canopy_crop_fixed(2,2)))
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,3)); % precip extinction coefficient
    else
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,2));
    end
    
    
    % PHOTOSYNTHESIS
    load './Temps/temp_variable.mat'...
        'fullpath_forcings' 'para_canopy_crop4' 'para_photosynthesisC3_crop4' 'para_photosynthesisC4_crop4' 'ph_type4'
    load (fullpath_forcings)
    
    PARAMS.Photosyn.ph_type(ii) = ph_type4;
    
    if PARAMS.Photosyn.ph_type(ii) == 1
        % C3 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC3_crop4(1,2)))
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop4(1,3));
        else
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop4(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop4(2,2)))
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop4(2,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop4(2,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop4(3,2)))
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop4(3,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop4(3,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop4(4,2)))
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop4(4,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop4(4,2))*ones(1,size(Ta_crop,1));
        end
        
        % C4 Photosynthesis Parameters
        PARAMS.Photosyn.Vmax_C4(ii) = NaN;
        PARAMS.Photosyn.Rd_C4(ii) = NaN;
        PARAMS.Photosyn.Q10_C4(ii) = NaN;
        PARAMS.Photosyn.kk_C4(ii) = NaN;
        PARAMS.Photosyn.theta_C4(ii) = NaN;
        PARAMS.Photosyn.beta_C4(ii) = NaN;
        PARAMS.Photosyn.al_C4(ii) = NaN;
    else
        % C3 Photosynthesis Parameters
        PARAMS.Photosyn.beta_ph_C3(ii) = NaN;
        PARAMS.Photosyn.Vcmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Jmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Rd25{ii} = NaN;
        
        % C4 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC4_crop4(1,2)))
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop4(1,3));
        else
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop4(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop4(2,2)))
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop4(2,3));
        else
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop4(2,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop4(3,2)))
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop4(3,3));
        else
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop4(3,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop4(4,2)))
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop4(4,3));
        else
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop4(4,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop4(5,2)))
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop4(5,3));
        else
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop4(5,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop4(6,2)))
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop4(6,3));
        else
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop4(6,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop4(7,2)))
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop4(7,3));
        else
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop4(7,2));
        end
    end
    
    PARAMS.Photosyn.Oi(ii) = 210;                                     % Intercellular oxygen concentration [mmol / mol]
    
    
    % STOMATAL CONDUCTANCE
    load './Temps/temp_variable.mat'...
        'para_conductance_crop4'
    
    % Ball-Berry
    if isempty(cell2mat(para_conductance_crop4(1,2)))
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop4(1,3)); % slope parameter in BB model [-] (Leakey: m=10.6 (ambient); m=10.9 (elevated))
    else
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop4(1,2));
    end
    if isempty(cell2mat(para_conductance_crop4(2,2)))
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop4(2,3)); % intercept parameter in BB model [mol/m^2/s] (Leakey: b=0.008 (ambient); b=0.007 (elevated))
    else
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop4(2,2));
    end
    
    % (Tuzet et al, PCE 2003)
    if isempty(cell2mat(para_conductance_crop4(3,2)))
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop4(3,3)); % sensitivity parameter for initial decrease in leaf potential [-]
    else
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop4(3,2));
    end
    if isempty(cell2mat(para_conductance_crop4(4,2)))
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop4(4,3)); % leaf potential at which half of the hydraulic conductance is lost [MPa]
    else
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop4(4,2));
    end
    
    % Conductivities
    if isempty(cell2mat(para_conductance_crop4(5,2)))
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop4(5,3)); % radial conductivity of the root system [s^-1]
    else
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop4(5,2));
    end
    if isempty(cell2mat(para_conductance_crop4(6,2)))
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop4(6,3)); % axial specific conductivity of the root system [mm/s]
    else
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop4(6,2));
    end
    if isempty(cell2mat(para_conductance_crop4(7,2)))
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop4(7,3));
    else
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop4(7,2));
    end
    
    if isempty(cell2mat(para_conductance_crop4(8,2)))
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop4(8,3));
    else
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop4(8,2));
    end
    
    
    % ROOT SOIL
    % HR
    load './Temps/temp_variable.mat'...
        'HR'
    SWITCHES.HR_on(ii) = HR;                                          % Hydraulic Redistribution 1. Yes, 0. No
    
    % root structure
    load './Temps/temp_variable.mat'...
        'set_para_root4'
    PARAMS.Soil.z50f(ii) = cell2mat(set_para_root4(2,2)); % From Ameriflux site description
    PARAMS.Soil.z95f(ii) = cell2mat(set_para_root4(3,2)); % From Ameriflux site description
    PARAMS.Soil.z50t(ii) = cell2mat(set_para_root4(2,2));
    PARAMS.Soil.z95t(ii) = cell2mat(set_para_root4(3,2));
    PARAMS.Soil.maxrootdepth(ii) = cell2mat(set_para_root4(1,2));
    
end

% Dongkook - LOAD PARAMETERS_species_bondville_soy End




nl_can = PARAMS.CanStruc.nl_can;
nspecies=PARAMS.CanStruc.nspecies;

%*************************************************************************%
%**************************** Flux / Met Data ****************************%
%******************************************************* *****************%

% % LOAD DATA
% if bon_type == 5
%     shdry = 0;
%     shring = 1;
%     if SWITCHES.CLIMATECHANGE == 1
%         shelev = 1;
%     else
%         shelev = 0;
%     end
%     [ELEV, LAT, LONG, ...
%         decdoy, year, doy, hour, month, ZEN_in, ...
%         U_in, Rg_in, ustar_in, PPT_in, Ta_in, ...
%         ea_in, VPD_in, Pa_in, HR_in,...
%         LAI1_in, LAI_in] = ...
%         LOAD_DATA_BONDVILLE_DRYFACE(years, doys, shdry, shelev, shring);
%     load('./Data_Files/BONDVILLE/SOYFACE/IC_sim_drysoy');
%     GPP1_in = nan;
%     NEE1_in = nan;
%     LE_in = nan;
%     H_in = nan;
% elseif bon_type == 6
%     shdry = 0;
%     shring = 1;
%     if SWITCHES.CLIMATECHANGE == 1
%         shelev = 1;
%     else
%         shelev = 0;
%     end
%     [ELEV, LAT, LONG, ...
%         decdoy, year, doy, hour, month, ZEN_in, ...
%         U_in, Rg_in, ustar_in, PPT_in, Ta_in, ...
%         ea_in, VPD_in, Pa_in, HR_in,...
%         LAI1_in, LAI_in] = ...
%         LOAD_DATA_BONDVILLE_DRYFACE_hh(years, doys, shdry, shelev, shring);
%     load('./Data_Files/BONDVILLE/SOYFACE/IC_sim_drysoy');
%     GPP1_in = nan;
%     NEE1_in = nan;
%     LE_in = nan;
%     H_in = nan;
% else
%     [ELEV,  LAT,        LONG,       decyear,        decdoy,     year,...
%         doy,        hour,       ZEN_in,         LAI_in,     Rg_in,...
%         Ta_in,      VPD_in,     PPT_in,         U_in,       ustar_in,...
%         Pa_in,      ea_in,      Fc_in,          LE_in,      H_in,...
%         G_in,       Fc_qc,      LE_qc,          H_qc,       Tskin_in,...
%         Ts4_in,     Ts8_in,     Ts16_in,        Ts32_in,    Ts64_in,...
%         Ts128_in,   SWC10_in,   SWC20_in,       SWC30_in,   SWC40_in,...
%         SWC50_in,   SWC60_in,   SWC100_in,      Rgup_in,    LWdn_in,...
%         LWup_in,    Rn_in] = ...
%         LOAD_DATA_SoyAmeriflux(years, doys, PARAMS, CONSTANTS, datafile, bon_type);
%     load('./Data_Files/BONDVILLE/SOYFACE/IC_sim_drysoy');
%
% end



% Dongkook - LOAD_DATA_SoyAmeriflux


% % Load Soy (Bondville) Ameriflux data
% %   - gapfilled .mat data file name is hardcoded
%
%   load(datafile);
%     laimin = 3;
%     inds = find(ismember(year_bon, years) & ...
%                 ismember(doy_bon, doys)& ...
%                 LAI_bon >= laimin);
%
%
%     decyear = decyear_bon(inds);
%     decdoy = decdoy_bon(inds);
%     year = year_bon(inds);
%     doy = doy_bon(inds);
%     hour = hour_bon(inds);
%     ZEN = ZEN_bon(inds);
%     LAI = LAI_bon(inds);
%
%     Ta = Ta_bon(inds);
%     VPD = VPD_bon(inds);
%     PPT = PPT_bon(inds);
%     U = U_bon(inds);
%     ustar = ustar_bon(inds);
%     Pa = Pa_bon(inds);
%     ea = ea_bon(inds);
%
%     Fc = Fc_bon(inds);
%     LE = LE_bon(inds);
%     H = H_bon(inds);
%     Hg = Hg_bon(inds);
%
%     Fc_qc = Fc_qc_bon(inds);
%     LE_qc = LE_qc_bon(inds);
%     H_qc = H_qc_bon(inds);
%
%     Tskin = Tskin_bon(inds);
%     Ts4 = Ts4_bon(inds);
%     Ts8 = Ts8_bon(inds);
%     Ts16 = Ts16_bon(inds);
%     Ts32 = Ts32_bon(inds);
%     Ts64 = Ts64_bon(inds);
%     Ts128 = Ts128_bon(inds);
%
%     SWC10 = SWC10_bon(inds);
%     SWC20 = SWC20_bon(inds);
%     SWC30 = SWC30_bon(inds);
%     SWC40 = SWC40_bon(inds);
%     SWC50 = SWC50_bon(inds);
%     SWC60 = SWC60_bon(inds);
%     SWC100 = SWC100_bon(inds);
%
%     Rg = Rg_bon(inds);
%     Rgout = Rgout_bon(inds);
%     LWin = LWin_bon(inds);
%     LWout = LWout_bon(inds);
%     Rn = Rn_bon(inds);
%
% % Calculate Vapor Variables
%     aa = 0.611;  % [kPa]
%     bb = 17.502;
%     cc = 240.97; % [C]
%     esat = aa*exp((bb*Ta)./(cc+Ta));
%     ean = esat - VPD;
%     hr = ea./esat;
%
% % Data Corrections
%     uinds = find(U<0.1);
%     U(uinds) = 0.1;
%
% % Calculate ustar for missing periods
%     binds = find(isnan(ustar) | ustar<=0);
%     vonk = 0.41;
%     hobs = 1.5;
%     z0 = 0.005;
%     hcan =1;
%     d0 = 2/3 * hcan;
%     ustar(binds) = vonk.*U(binds)./(log((hobs-d0)/z0));
%
% % Adjust U from measurement height to canopy height (Brutsaert, 4.2)
%     U = U - (ustar/vonk).*log(hobs/hcan);
%     uinds = find(U<0.1);
%     U(uinds) = 0.1;
%
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% % FILL Measured Hg with Linear fit of Hg to Rg
% % --> Using Hg to force soil heat transfer model
%     ginds = find(~isnan(Hg) & ~isnan(Rg));
%     [lps] = polyfit(Rg(ginds), Hg(ginds), 1);
%     aa = lps(1); bb = lps(2);
%     binds = find(isnan(Hg));
%     Hg(binds) = aa*Rg(binds) + bb;
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% LOAD DATA
decyear = nan;

load './Temps/temp_variable.mat'...
    'fullpath_forcings' 'lat_face' 'long_face' 'elev_face'
load (fullpath_forcings)

% Organize the data
inds = find(ismember(year_crop, Each_year) & ismember(floor(round((doy_crop+1).*10^10)./10^10), doys));

ELEV      = elev_face;
LAT       = lat_face;
LONG      = long_face;
decyear   = year_crop(inds);
decdoy    = doy_crop(inds);

year      = year_crop(inds);
doy       = doy_crop(inds);
hour      = hour_crop(inds);
ZEN_in    = ZEN_crop(inds);
Rg_in     = Rg_crop(inds);

Ta_in     = Ta_crop(inds);
VPD_in    = VPD_crop(inds);
PPT_in    = PPT_crop(inds);
U_in      = U_crop(inds);
ustar_in  = ustar_crop(inds);

Pa_in     = Pa_crop(inds);
ea_in     = ea_crop(inds);

% Dongkook: Calucalated in model_forcing
% % Calculate Vapor Variables
%     aa = 0.611;  % [kPa]
%     bb = 17.502;
%     cc = 240.97; % [C]
%     esat = aa*exp((bb*Ta)./(cc+Ta));
%     ean = esat - VPD;
%     hr = ea./esat;
%

% Data Corrections
uinds = find(U_in<0.1);
U_in(uinds) = 0.1;

% Calculate ustar for missing periods
binds = find(isnan(ustar_in) | ustar_in<=0);
vonk = CONSTANTS.vonk;
hobs = PARAMS.CanStruc.hobs;
z0 = PARAMS.CanStruc.z0;
hcan =PARAMS.CanStruc.hcan;
d0 = PARAMS.CanStruc.d0;
ustar_in(binds) = vonk.*U_in(binds)./(log((hobs-d0)/z0));

% Adjust U from measurement height to canopy height (Brutsaert, 4.2)
U_in = U_in - (ustar_in/vonk).*log(hobs/hcan);
uinds = find(U_in<0.1);
U_in(uinds) = 0.1;

% Dongkook: No idea where it is used
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% % FILL Measured Hg with Linear fit of Hg to Rg
% % --> Using Hg to force soil heat transfer model
%     ginds = find(~isnan(Hg) & ~isnan(Rg));
%     [lps] = polyfit(Rg(ginds), Hg(ginds), 1);
%     aa = lps(1); bb = lps(2);
%     binds = find(isnan(Hg));
%     Hg(binds) = aa*Rg(binds) + bb;
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% ustar_in  =
% Fc_in     =
% LE_in     =
% H_in      =
% LAI_in    =
% G_in      =
% Fc_qc     =
% LE_qc     =
% H_qc      =
% Tskin_in  =
% Ts4_in    =
% Ts8_in    =
% Ts16_in   =
% Ts32_in   =
% Ts64_in   =
% Ts128_in  =
% SWC10_in  =
% SWC20_in  =
% SWC30_in  =
% SWC40_in  =
% SWC50_in  =
% SWC60_in  =
% SWC100_in =
% Rgup_in   =
% LWdn_in   =
% LWup_in   =
% Rn_in     =

if (kspecies == 1)
    LAI_in(:,1)= LAI_species1_crop(inds);
elseif (kspecies == 2)
    LAI_in(:,1)= LAI_species1_crop(inds);
    LAI_in(:,2)= LAI_species2_crop(inds);
elseif (kspecies == 3)
    LAI_in(:,1)= LAI_species1_crop(inds);
    LAI_in(:,2)= LAI_species2_crop(inds);
    LAI_in(:,3)= LAI_species3_crop(inds);
elseif (kspecies == 4)
    LAI_in(:,1)= LAI_species1_crop(inds);
    LAI_in(:,2)= LAI_species2_crop(inds);
    LAI_in(:,3)= LAI_species3_crop(inds);
    LAI_in(:,4)= LAI_species4_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 1)
    LAI_in(:,1)= LAI_species1_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 2)
    LAI_in(:,1)= LAI_species2_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 3)
    LAI_in(:,1)= LAI_species3_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 4)
    LAI_in(:,1)= LAI_species4_crop(inds);
end

% Dongkook - LOAD_DATA_SoyAmeriflux End




% SET UP CANOPY GRID and LEAF AREA DISTRIBUTION
load './Temps/temp_variable.mat'...
    'dat_LAD1' 'dat_LAD2' 'dat_LAD3' 'dat_LAD4'

% [znc, zhc, dzc] = CANOPY_GRID(PARAMS);
% VERTSTRUC.znc = znc;
% VERTSTRUC.zhc = zhc;
% VERTSTRUC.dzc = dzc;

znc = (dat_LAD1(:,1))';                                               % Middle of each canopy grid [m]
zhc = (dat_LAD1(:,1)+(dat_LAD1(2,1)-dat_LAD1(1,1))/2)';               % Top of each canopy grid [m]
%dzc = zhc(1)./2;
dzc = zhc(1);

VERTSTRUC.znc = znc;
VERTSTRUC.zhc = zhc;
VERTSTRUC.dzc = dzc;



% if bon_type == 1
%     [LADnorm] = LAD_Corn_Boedhram(PARAMS, VERTSTRUC);
% elseif (bon_type == 2 || bon_type == 5 || bon_type == 6)
%     [LADnorm] = LAD_Dermody(VERTSTRUC);
% elseif bon_type == 3
%     [LADnorm] = LAD_Kromdijk(VERTSTRUC);
% elseif bon_type == 4
%     [LADnorm] = LAD_Madakadze(VERTSTRUC);
% end
% LADnorm_all = LADnorm;


if (kspecies == 1)
    LADnorm_all(:,1)= dat_LAD1(:,2);
elseif (kspecies == 2)
    LADnorm_all(:,1)= dat_LAD1(:,2);
    LADnorm_all(:,2)= dat_LAD2(:,2);
elseif (kspecies == 3)
    LADnorm_all(:,1)= dat_LAD1(:,2);
    LADnorm_all(:,2)= dat_LAD2(:,2);
    LADnorm_all(:,3)= dat_LAD3(:,2);
elseif (kspecies == 4)
    LADnorm_all(:,1)= dat_LAD1(:,2);
    LADnorm_all(:,2)= dat_LAD2(:,2);
    LADnorm_all(:,3)= dat_LAD3(:,2);
    LADnorm_all(:,4)= dat_LAD4(:,2);
elseif (kspecies == 5 && Sim_species_con == 1)
    LADnorm_all(:,1)= dat_LAD1(:,2);
elseif (kspecies == 5 && Sim_species_con == 2)
    LADnorm_all(:,1)= dat_LAD2(:,2);
elseif (kspecies == 5 && Sim_species_con == 3)
    LADnorm_all(:,1)= dat_LAD3(:,2);
elseif (kspecies == 5 && Sim_species_con == 4)
    LADnorm_all(:,1)= dat_LAD4(:,2);
end

% ALLOCATE MEMEORY FOR NVINDS
nvinds_all = cell(nspecies,1);
vinds_all = cell(nspecies,1);


% Dongkook Woo - Edit End