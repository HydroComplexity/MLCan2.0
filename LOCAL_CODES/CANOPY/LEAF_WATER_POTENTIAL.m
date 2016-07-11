function [psil_MPa] = LEAF_WATER_POTENTIAL (VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, cntspecies)


%=========================================================================
%   This function is used to calculates Leaf Water Potential at each      %  
%   layer using Ohm's law.                                                %
%   See Drewry et al, 2009, Part B Online Supplement, Eqns (13,14)        %
%
%   Written By: Darren Drewry
%   All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES       % VARIABLES structure
%       PARAMS          % PARAMS structure
%       VERTSTRUC       % VERTSTRUC structure
%       CONSTANTS       % CONSTANTS structure
%       cntspecies      % [] number of species
%------------------------- Output Variables ------------------------------
%       psil_MPa        % [MPa] Leaf Water Potential
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
%
    % VARIABLES
    TR_sun      = VARIABLES.CANOPY.TR_sun(:,cntspecies);                   % [mm/s/unit LAI] transpiration PER UNIT LEAF AREA sunlit                   
    TR_shade    = VARIABLES.CANOPY.TR_shade(:,cntspecies);                 % [mm/s/unit LAI] transpiration PER UNIT LEAF AREA shade     
    rpp_wgt     = VARIABLES.ROOT.rpp_wgt(:,cntspecies);                    % [mm] Root pressure potential weighted by root distribution

    % VERTSTRUC
    znc         = VERTSTRUC.znc;                                           % [m] Height of canopy levels 

    % PARAMS
    Rp          = PARAMS.StomCond.Rp(cntspecies);                          % [MPa / m / s] Plant resistance to water flow

    % CONSTANTS   
    grav        = CONSTANTS.grav;                                          % % [m / s^2] Gravity Acceleration
    dtime       = CONSTANTS.dtime;                                         % [s] Time Step 
    mmH2OtoMPa  = CONSTANTS.mmH2OtoMPa;                                    % [] Conversion Factor from mmH2O to MPa
%    
%*************************************************************************%
%% <<<<<<<<<<<<<<<<<<<<<<< END OF DE-REFERENCE BLOCK >>>>>>>>>>>>>>>>>>> %%
%*************************************************************************%
%% 
    rho_kg      = 1;                                                        % [kg / m^3]
%
    znc         = znc(:);                                                   % [m]
%
    TR          = (TR_sun + TR_shade);                                      % [W/LAI/s]
%
    TR_m        = TR / 1000;                                                % [m/s/unit LAI]
%
    rpp_wgt_MPa = rpp_wgt * mmH2OtoMPa;                                     % [MPa]
%
    psil_MPa    = rpp_wgt_MPa - TR*Rp - (rho_kg*grav*znc)./10^6;            % [MPa]
%
%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%% <<<<<<<<<<<<<<<<<<<<<<<<< END OF FUNCTION >>>>>>>>>>>>>>>>>>>>>>>>>>>>%%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%    
