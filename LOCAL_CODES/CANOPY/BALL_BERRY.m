function [gsv, Ci, Hs, Cs] = BALL_BERRY (VARIABLES, PARAMS, VERTSTRUC, sunlit,cntspecies, SWITCHES, FORCING)

%=========================================================================
%   Calculates stomatal conductance using the Ball-Berry equation
%
% Written By: Darren Drewry, Modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES        % VARIABLES structure
%       PARAMS           % PARAMS structure
%       VERTSTRUC        % PARAMS structure
%       sunlit           % sunlit fraction
%       cntspecies       % Number of species
%       SWITCHES         % SWITCHES structure
%       FORCING          % PARAMS structure
%------------------------- Output Variables ------------------------------
%       gsv              % [mol/m^2 la/s] gsv - stomatal conductance to vapor transport  
%       Ci               % [umol/mol] Internal leaf concentration of CO2   
%       Hs               % [] Relative Humidity in the leaf boundary layer boundary layer   
%       Cs               % [umol/mol] Concentration of CO2 in leaf boundary layer                
%========================================================================              
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES FOR EACH SPECIES 
    %**********************************************************************
    if (sunlit)
        An =  VARIABLES.CANOPY.An_sun(:,cntspecies);                       % [W/m^2 leaf area] Rabs - Total absorbed radiation Sunlit Fraction     
        Tl_C =  VARIABLES.CANOPY.Tl_sun(:,cntspecies);                     % [C] Leaf temperature canopy sunlit fraction 
        gsv =  VARIABLES.CANOPY.gsv_sun(:,cntspecies);                     % [mol/m^2 la/s] gsv - stomatal conductance to vapor transport 
        fsvg =  VARIABLES.CANOPY.fsvg_sun(:,cntspecies);                   % [mol/m^2/s] Stomatal Conductance for sunlit fraction for all layers and species 
        fsvm =  VARIABLES.CANOPY.fsvm_sun(:,cntspecies);                   % [] m Tuzet factor to Ball Berry Model for sunlit fraction, for all layers and species 
        gbv = VARIABLES.CANOPY.gbv_sun(:,cntspecies);                      % [mol/m^2/s] Boundary Layer conductance to vapor in sunlit, for all layers and species 
    else
        An =  VARIABLES.CANOPY.An_shade(:,cntspecies);                     % [W/m^2 leaf area] Rabs - Total absorbed radiation Sunlit Fraction  
        Tl_C =  VARIABLES.CANOPY.Tl_shade(:,cntspecies);                   % [C] Leaf temperature canopy shade fraction 
        gsv =  VARIABLES.CANOPY.gsv_shade(:,cntspecies);                   % [mol/m^2 la/s] gsv - stomatal conductance to vapor transport 
        fsvg =  VARIABLES.CANOPY.fsvg_shade(:,cntspecies);                 % [mol/m^2/s] Stomatal Conductance for sunlit fraction for all layers and species 
        fsvm =  VARIABLES.CANOPY.fsvm_shade(:,cntspecies);                 % [] m Tuzet factor to Ball Berry Model for sunlit fraction, for all layers and species 
        gbv =  VARIABLES.CANOPY.gbv_shade(:,cntspecies);                   % [mol/m^2/s] Boundary Layer conductance to vapor in sunlit, for all layers and species 
    end
        
    Ciprev =  VARIABLES.CANOPY.Ci_sun(:,cntspecies);                       % [umol/mol] Internal leaf concentration of CO2 in previous time step

    mslope = PARAMS.StomCond.mslope(cntspecies);                           % [] Ball Berry parameter. slope. 
    bint = PARAMS.StomCond.bint(cntspecies);                               % [] Ball Berry parameter. Intercept. 
    %*********************************************************************       
    nvinds_all = VERTSTRUC.nvinds_all;                                     % [] Indicator of Layers where LAI is zero of lower for all species 
    nvinds = nvinds_all{cntspecies};                                       % [] Indicator of Layers where LAI is zero of lower 
    Ca = VARIABLES.CANOPY.CAz;                                             % [umol/mol] Atmospherit concentration of CO2  
    ea = VARIABLES.CANOPY.EAz;                                             % [kPa] Vapor Pressure Air, for all layers     
    
%*************************************************************************


% Fick's Law (18c - Sellers et al, 1992)
    Cs = Ca - (1.37 * An) ./ gbv;    
% Ball / Berry
    estarTl = 0.611 * exp( (17.502*(Tl_C)) ./ (Tl_C + 240.97) ); % [kPa]
    ei = estarTl; % Assume inter-cellular space saturated at Tl            
    es = (gsv.*ei + gbv.*ea) ./ (gsv + gbv);
    Hs = es ./ estarTl; 
    
    gsv = fsvm.*(mslope.*An.*Hs./Cs) + (bint);

% Internal CO2 Concentration
    Ci = Ca - (1.37 * An)./gbv - (1.6 * An)./gsv;
%    ind=gsv==0;
%    Ci(ind)=Ciprev(ind);
% No net uptake
    binds = find(An<=0);
    Ci(binds) = Ca(binds);
    gsv(binds) = (bint);
%*************************    
    indg = gsv<bint; 
    gsv(indg) = bint;  
%*************************    
% Non-vegetated levels    
    gsv(nvinds) = 0;
    Ci(nvinds) = 0;
        
    gsv = gsv(:);
    Ci = Ci(:);