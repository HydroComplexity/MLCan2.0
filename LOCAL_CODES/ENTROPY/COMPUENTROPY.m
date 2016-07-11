
function [SSresults] = ...
    COMPUENTROPY (SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,... 
                  SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
                  SWout, fdiff,LWabs_canM, LWabs_soilM, LWemit_soil, LWemit_sun, LWemit_shade,...
                         LWout, Tsurf, FORCING,SWITCHES, CONSTANTS, PARAMS, VARIABLES, VERTSTRUC)

%=========================================================================
% This code computes the fluxes of entropy in the ecosystem. 
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       SWcandir_in                    % [W/m2] Direct Incoming Shortwave Radiation in canopy
%       SWcandir_out                   % [W/m2] Direct Outgoing Shortwave Radiation in canopy
%       SWcandif_in                    % [W/m2] Diffuse Incoming Shortwave Radiation in canopy 
%       SWcandif_out                   % [W/m2] Diffuse Outgoing Shortwave Radiation in canopy
%       SWsoildir_in                   % [W/m2] Direct Incoming Shortwave Radiation in soil
%       SWsoildir_out                  % [W/m2] Direct Outgoing Shortwave Radiation in soil
%       SWsoildif_in                   % [W/m2] Diffuse Incoming Shortwave Radiation in soil
%       SWsoildif_out                  % [W/m2] Diffuse Outgoing Shortwave Radiation in soil
%       SWout                          % [W/m2] Total Outgoing Shortwave Radiation 
%       fdiff                          % [] Canopy Fraction that receives diffuse
%       LWabs_canM                     % [W/m2] Matrix of Absorptin of LW radiation in canopy, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWabs_soilM                    % [W/m2] Vector of Absorptin of LW radiation in the soil, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWemit_soil                    % [W/m2] LW radiation emitted from soil
%       LWemit_sun                     % [W/m2] LW radiation emitted from canopy sunlit
%       LWemit_shade                   % [W/m2] LW radiation emitted from canopy sunlit
%       LWout                          % [W/m2] Total LW radiation out 
%       Tsurf                          % [C] Soil Surface Temperature
%       FORCING                        % FORCING structure
%       SWITCHES                       % SWITCHES structure
%       CONSTANTS                      % CONSTANTS structure
%       PARAMS                         % PARAMS structure
%       VARIABLES                      % VARIABLES structure
%       VERTSTRUC                      % VERTSTRUC structure
%------------------------- Output Variables ------------------------------
%       SSresults                      % SSresults structure 
% 
%========================================================================
                              
% De reference block
%*********************************************************************
% CONSTANTS
boltz = CONSTANTS.boltz;         % Stefan-Boltzmann constant [W/m^2/K^4]

% SWITCHES
entropymethod = SWITCHES.entropymethod;

% FORCING
Rg=FORCING.Rg;
Tatop = FORCING.Ta;        
LWin = FORCING.LWdn;

% VERTSTRUC
fLAIz = VERTSTRUC.fLAIz;

% PARAMS
nspecies = PARAMS.CanStruc.nspecies;

% VARIABLES
zicesl = VARIABLES.SOIL.zicesl;

%**************************************************************************
% COMPUTATION OF ENTROPY DUE TO SHORTWAVE 
            
            
 [SSresults] = ENTROPY_SW (SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
                                          SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
                                          SWout, fdiff, entropymethod, Rg, zicesl, PARAMS);
                                               
%**************************************************************************                                   
% COMPUTATION OF ENTROPY DUE TO LONGWAVE 


[SSresults] = ENTROPY_LW (LWabs_canM, LWabs_soilM, LWemit_soil, LWemit_sun, LWemit_shade,...
                                   LWin, LWout, Tatop, Tsurf, boltz, entropymethod, zicesl, ...
                                   PARAMS, VARIABLES, CONSTANTS, SWITCHES, SSresults);


%**************************************************************************                                   
% COMPUTATION OF ENTROPY DUE TO 	HEAT FLUXES 
                               
[SSresults] = ENTROPY_HEAT (VARIABLES, SWITCHES, SSresults);                           


%**************************************************************************                                   

[SSresults] = ENTROPY_PHO (VARIABLES, SWITCHES, VERTSTRUC, PARAMS, SSresults);

[SSresults] = ENTROPY_results(SSresults, VERTSTRUC, PARAMS, VARIABLES);

[SSresults] = ENTROPY_H2OINF (VARIABLES, VERTSTRUC, CONSTANTS, PARAMS, SSresults);                               

