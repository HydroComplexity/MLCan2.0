function [fsvg fsvm] = Tuzet_Function (VARIABLES, PARAMS, sunlit, cntspecies)

%=========================================================================
%   Calculates the stomatal conductance reduction
%   function due to plant hydraulic constraint from Tuzet et al (PCE, 2003)
%
% Written By: Darren Drewry
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING         % FORCING structure
%       PARAMS          % PARAMS structure
%       sunlit          % sunlit indicator 
%       cntspecies      % [] number of species
%------------------------- Output Variables ------------------------------
%       fsvg            %[] Applied factor for the intercept part in the Ball Berry Model
%       fsvm            %[] Applied factor for the slope part in the Ball Berry Model
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************

    % VARIABLES FOR EACH SPECIES 
    %**********************************************************************
    if (sunlit)
        psil_MPa = VARIABLES.CANOPY.psil_sun(:,cntspecies);                % [MPa] Leaf Water Potential
    else
        psil_MPa = VARIABLES.CANOPY.psil_shade(:,cntspecies);              % [MPa] Leaf Water Potential     
    end
    
    % PARAMS
    nlayers = PARAMS.CanStruc.nl_can;                                      % [] Number of layers in canopy 
    sfm = PARAMS.StomCond.sfm(cntspecies);                                 % Sensitivity parameter for initial decrease in leaf potential [-]  
    psifm = PARAMS.StomCond.psifm(cntspecies);                             % Leaf potential at which half of the hydraulic conductance is lost [MPa] 
%**********************************************************************
%**************************************************************************
    
%     fsvg = (1 + exp(sfg*psifg)) ./ (1 + exp(sfg*(psifg - psil_MPa)));
%     if cntspecies == 1
%         if abs(psil_MPa) > 2
%              fsvm = (1 + exp(sfm*psifm)) ./ (1 + exp(sfm*(psifm - psil_MPa)));
%         else
%              fsvm=ones(length(psil_MPa),1);
%         end
%            fsvg = min(0.1,fsvg);
%     else
        fsvm = (1 + exp(sfm*psifm)) ./ (1 + exp(sfm*(psifm - psil_MPa)));
        fsvg = ones(nlayers,1);
%    end
    