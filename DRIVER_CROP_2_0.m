function DRIVER_CROP_2_0
% Canopy-Root-Soil-Atmosphere Exchange Model
%
%   Written By : Darren Drewry (dtd2@illinois.edu)
%              : Juan Quijano (quijano2@illinois.edu)
%              : Dongkook Woo (dwoo5@illinois.edu)

% STRUCTURES:
%   SWITCHES --> model conditional switches
%   VERTSTRUC --> variables describing vertical structure of canopy & soil
%   VARIABLES.CANOPY --> canopy variables
%            .SOIL --> soil variables
%            .ROOT --> root variables
%   FORCING --> holds current timestep forcing variables
%   CONSTANTS --> site independent constants, unit conversions
%   PARAMS
%         .CanStruc --> canopy structural parameters
%         .Rad --> radiation parameters
%         .Photosyn --> photosynthesis paramters
%         .StomCond --> stomtatal conductance parameters
%         .Resp --> ecosystem respiration parameters
%         .MicroEnv --> canopy microenvironment parameters
%         .Soil --> soil paramters
%
%**************************************************************************
%                           USER SPECIFICATIONS
%**************************************************************************
load './Temps/temp_variable.mat' ...
    'DOY_start' 'DOY_end'
doys = [DOY_start:DOY_end];

%******************************  SWITCHES   *******************************
load './Temps/temp_variable.mat'...
    'Turbulence' 'set_para_root1'  'vanGen' 'RHC' 'Soil_heat' ...
    'Temp_Elev' 'Temp_Elev_con' 'CO2_Ambient' 'CO2_Elev'  ...
    'CO2_Elev_con' 'Soil_nutrient' 'litter_depth'

SWITCHES.plants = 1;              %  Allow the presence or no of plants. 1. With plants 2. No plants (bare soil, no roots)
SWITCHES.LT = 0;                  %  Run long term dynamics with stochastically generated data 1. Yes 0. No
SWITCHES.turb_on = Turbulence;    %  1 = scalar profiles resolved, otherwise not resolved
if litter_depth > 0
    SWITCHES.litter = 1;          %  1 = Include litter dynamics. 0 does not include litter dynamics
else
    SWITCHES.litter = 0;          %  1 = Include litter dynamics. 0 does not include litter dynamics
end
SWITCHES.littertype = 1;          %  [1 LEs parameterized with aerodynamic resitance

%  2 LEs parameterized with vapor litter diffusivity]
SWITCHES.soilevap = 1;            %  1 = Include soil evaporation. 0 does not include soil evaporation
if set_para_root1{4,2} > 0
    SWITCHES.cutroots = 1;        %  0 if no cut at all
                                  %  (nlc) 1 if cut all the time same layers
                                  %  2 if cut with a threshold
                                  %  3 if cut with a embolism curve (PLC)
                                  %  4 combination of 1, 2 and 3
elseif set_para_root1{4,2} == 0
    SWITCHES.cutroots = 0;        %  0 if no cut at all
                                  %  (nlc) 1 if cut all the time same layers
                                  %  2 if cut with a threshold
                                  %  3 if cut with a embolism curve (PLC)
                                  %  4 combination of 1, 2 and 3
end

SWITCHES.rtcond  = 0;             %  0. Juan's Method, 1. Amenu Method
SWITCHES.rhc = RHC;               %  1 = linearly increasing root hydraulic conductivity with depth
SWITCHES.ns=1;                    %  Numerical Scheme. 1 = Use Implicit, 0 = Use explicit

SWITCHES.soilheat_on = Soil_heat; %  Compute the heat equation 1 = Yes, 0 = No
SWITCHES.canstorheat = 1;         %  Include storage of heat in leaves 1 = Yes, 0 = No
SWITCHES.useG_on=0;               %  Use of G measured in the soil energy balance instead of compute it
SWITCHES.useTs_on=0;              %  Use of Ts (soil temperature) [soil energy balance]

SWITCHES.entropy_on = 1;          %  Compute Entropy, 1 = Yes, 0 = No.
SWITCHES.entropymethod = 2;       %  Compute Entropy Method for LW entropy. 1. Landsberg 1979. 2. Smith 2001. 3. Clausius

SWITCHES.save_on = 1;             %  1 = save stored variables to .mat file, otherwise no save performed
SWITCHES.plots_on = 0;

SWITCHES.fsv_off = 0;             %  1 = hydraulic constraint turned OFF

SWITCHES.temp_change = Temp_Elev; %  0 = +0 C, 1 = 1 C,  2 = 2 C

% All of the CN model switch and parameters are in core_N
SWITCHES.soilCN_on = Soil_nutrient;%  Compute CN dynamics. 1 = Yes, 0 = No.
SWITCHES.Pedofunctions = 0;       %  1. Use pedo transfer functions with organic matter 0. Do not use
SWITCHES.vanGen = vanGen;         %  1. vanGenuchten, 0. Brooks Corey / Soil moisture chracteristic.

%**************************************************************************
% Code Library Paths
addpath('./LOCAL_CODES/CANOPY/');
addpath('./LOCAL_CODES/ENTROPY/');
addpath('./LOCAL_CODES/ROOT_SOIL/');
addpath('./LOCAL_CODES/ROOT_SOIL/IMPLICIT');
addpath('./LOCAL_CODES/ROOT_SOIL/CN_MODEL');
addpath('./LOCAL_CODES/NUMERICAL/');
addpath('./LOCAL_CODES/NUMERICAL/OTHERS');
addpath('./LOCAL_CODES/PLOTTING/');

%************************** LOAD INFORMATION  *****************************
year = nan;                       % Initialize year, otherwise it is recognized as function
LOAD_SITE_INFO;

ybeginds = nan(length(Each_year),1);
yendinds = nan(length(Each_year),1);

for yy = 1:length(Each_year)
    ybeginds(yy) = find(year==Each_year(yy), 1, 'first');
    yendinds(yy) = find(year==Each_year(yy), 1, 'last');
end
%************************** CLIMATE CHANGE ********************************
if (SWITCHES.temp_change == 0)    % Ambient Temperature
    temp_change = 0;
elseif SWITCHES.temp_change == 1  % Elevated Temperature
    temp_change = Temp_Elev_con;
    Ta_in = Ta_in + temp_change;
end

CO2base_elevated = CO2_Elev_con;
CO2base_ambient = str2num(CO2_Ambient);
if (CO2_Elev == 1)                % Elevated CO2
    CO2base = CO2base_elevated;
else                              % Ambient CO2
    CO2base = CO2base_ambient;
end

%************************ CANOPY - SOIL - FUNCTION ************************
CANOPY_SOIL_COUPLER;
%**************************************************************************

%******************************* SAVE FUNCTION ****************************
load './Temps/temp_variable.mat'...
    'num_species'
if isstr(num_species) == 1
    Number_S = num_species;
else
    Number_S = num2str(num_species);
end
current_time=clock;
current_year=num2str(current_time(1));
current_month=num2str(current_time(2));
current_day=num2str(current_time(3));
current_hour=num2str(current_time(4));
current_min=num2str(current_time(5));

namefileout = ['./Result' '/Result_', Number_S, 'species', '_Saved_on_' ...
    , current_month, '.', current_day, '.', current_year, '.mat'];

% Make dianuual averages
MAKE_AVERAGES;

% Show message
if (length(Each_year)==1)
    timevect  = decdoy;
    timelabel = ['DOY: ', num2str(Each_year)];
else
    timevect  = [1:length(decyear)];
    %       timevect  = decyear;
    timelabel = ['Decimal Year'];
end

msgbox(['Simulation is done.', namefileout(10:end), ' is saved in the Result folder'],'MLCan Simulation');
save (namefileout)
%**************************************************************************


