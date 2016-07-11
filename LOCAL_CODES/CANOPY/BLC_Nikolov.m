function [gbv, gbh, gbv_forced, gbv_free, cnt] = BLC_Nikolov(VARIABLES, PARAMS, sunlit, cntspecies)

%=========================================================================
%   Computes Leaf Boundary Layer Conductance
%
%   Free and forced convective regimes are both computed, with free
%   convection becoming relevant under low wind conditions.
%
%   These equations are taken from (Nikolov, Massman, Schoettle),
%   Ecological Modelling, 80 (1995), 205-235
%   Written By: Darren Drewry, Modified by Juan Quijano
%   All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES        % VARIABLES structure
%       PARAMS           % PARAMS structure
%       sunlit           % sunlit fraction
%       cntspecies       % Number of species
%------------------------- Output Variables ------------------------------
%       gvb              % [mol/m^2/s] boundary layer conductance to vapor           
%       gvh              % [mol/m^2/s] boundary layer conductance to heat           
%       gvb_forced       % [mol/m^2/s] gvb forced           
%       gvb_free         % [mol/m^2/s] gvb free           
%========================================================================              
%*************************************************************************
%                          DE-REFERENCE BLOCK
%**************************************************************************
    % VARIABLES    
    if (sunlit)  
        Tl_in =  VARIABLES.CANOPY.Tl_sun(:,cntspecies);                    %[] Leaf Temperature sunlit fraction    
        gsv_in =  VARIABLES.CANOPY.gsv_sun(:,cntspecies);                  %[mol/m^2/s] Stomatal conductance for vapor sunlit fraction                
    else
        Tl_in =  VARIABLES.CANOPY.Tl_shade(:,cntspecies);                  %[] Leaf Temperature sunlit fraction 
        gsv_in =  VARIABLES.CANOPY.gsv_shade(:,cntspecies);                %[mol/m^2/s] Stomatal conductance for vapor sunlit fraction 
    end
    Ta_in = VARIABLES.CANOPY.TAz;                                          %[C] Air Temperature  
    Pa_in = VARIABLES.CANOPY.PAz;                                          %[kPa] Air pressure 
    ea_in = VARIABLES.CANOPY.EAz;                                          %[kpa] Vapor water potential i air 
    U = VARIABLES.CANOPY.Uz;                                               %[m/s] horizontal wind speed      

    % PARAMS
    ld = PARAMS.CanStruc.ld(cntspecies);                                   %[m] characteristic leaf dimension parameter
%          = leaf width for broadleaved vegetation
%          = needle diameter for conifers    
    lw = PARAMS.CanStruc.lw(cntspecies);                                   %[m] shoot diameter for conifers  
%          lw= 0 for broadleaved vegetation    
    
    leaftype = PARAMS.CanStruc.leaftype(cntspecies);                       %[] Type of leaf 1. Broadlife, 2. Needle  

%**************************************************************************
%**************************************************************************
% UNIT CONVERSIONS
    gsv = gsv_in / 41.4;    % [m/s]
    Tak = Ta_in + 273.15;   % [K]
    Tlk = Tl_in + 273.15;   % [K]
    Pa = Pa_in * 1000;      % [Pa]
    ea = ea_in * 1000;      % [Pa]
    
% SATURATION VAPOR PRESSURE at Tl    
    esTl = 1000 * 0.611*exp( (17.502*Tl_in) ./ (Tl_in+240.97) );    % [Pa]

% FORCED CONVECTION
    if (leaftype == 1)      % broadleaf
        cf = 1.6361 * 10^-3; %4.322 * 10^-3;
    elseif (leaftype == 2)  % evergreen needles
        cf = 0.8669 * 10^-3;%1.2035 * 10^-3;
    end    
        
    % [m/s]  Eqn (29) 
    gbv_forced = cf * Tak.^(0.56) .* ((Tak+120).*(U./ld./Pa)).^(0.5);
    

% FREE CONVECTION
    if (leaftype == 1)          % broadleaf
        ce = 1.6361 * 10^-3;
    elseif (leaftype == 2)      % evergreen needles
        ce = 0.8669 * 10^-3;
    end
    
    
% FORCED CONVECTION     
    gbv_free = gbv_forced;
    eb = (gsv.*esTl + gbv_free.*ea)./(gsv + gbv_free);  % Eqn 35

    Tvdiff = (Tlk ./ (1-0.378*eb./Pa)) - (Tak ./ (1-0.378*ea./Pa)); % Eqn 34

    gbv_free = ce * Tlk.^(0.56) .* ((Tlk+120)./Pa).^(0.5) .* (abs(Tvdiff)/lw).^(0.25); % Eqn 33
        
    
gbv_forced = gbv_forced * 41.4; % [mol/m^2/s]
gbv_free = gbv_free * 41.4;     % [mol/m^2/s]

gbv = max(gbv_forced, gbv_free);
gbh = 0.924 * gbv;

