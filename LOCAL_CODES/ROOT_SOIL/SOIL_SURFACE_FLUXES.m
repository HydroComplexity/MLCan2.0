function [Hs, LEs, Gs, RH, Ts, remain, VARIABLES] = SOIL_SURFACE_FLUXES(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES,FORCING)


%=========================================================================
% Solve surface energy balance based on formulation of (Hinzman et al, JGR 1998)
%
% Written by Darren Drewry, Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES       % VARIABLES structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%       CONSTANTS       % CONSTANTS structure
%       SWITCHES        % SWITCHES structure
%       FORCING         % CONSTANTS structure
%------------------------- Output Variables ------------------------------
%       Hs              % [W/m2] Sensible heat from Snow-Litter pack   
%       LEs             % [W/m2] Latent heat from Snow-Litter pack   
%       Gs              % [W/m2] Ground heat flux   
%       RH              % [] Relative Humidity of Top Soil Layer  
%       Ts              % [C] Temperature in the surface (Soil-Litter Interface)
%       Remain          % [W/m2] Error in soil energy balance
%       VARIABLES       % VARIABLES structure
% 
%========================================================================


%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES
    Ta1 = VARIABLES.CANOPY.TAz(1);                                          % [C] Atmospheric Temperature first layer
    ea1 = VARIABLES.CANOPY.EAz(1);                                          % [kPa] Vapor pressure at surface
    pa1 = VARIABLES.CANOPY.PAz(1);                                          % [kpa] Air pressure at surface 
    U1 = VARIABLES.CANOPY.Uz(1);                                            % [m/s] wind speed at surface       
    Rabs = VARIABLES.SOIL.Totabs_soil;                                      % [W / m2 / s] Radiation Energy Absorbed Soil
    volliq1 = VARIABLES.SOIL.volliq(1);                                     % [] Volumetric water content
    psis1 = VARIABLES.SOIL.smp(1);                                          % [MPa] Soil water potential of top layer
    
    % VERTSTRUC
    z1 = VERTSTRUC.znc(1);                                                  % [m] Height of bottom atmosphere node 
    dzs1 = VERTSTRUC.dzs(1);                                                % [m] Thickness of top soil layer
    TC1 = VERTSTRUC.TK_sol(1);                                              % [W / m / K] Thermal Conductivity firs soil layer
    porsl1 = VERTSTRUC.porsl(1);                                            % [] Porosity Top Soil Layer  
    
    %PARAMS
    epss = PARAMS.Rad.epss;                                                 % epss Soil emissivity
    z0 = PARAMS.Soil.z0;                                                    % [m] surface roughness length 
    
    %CONSTANTS
    vonk = CONSTANTS.vonk;                                                  % [] von Karman constant   
    Lv = CONSTANTS.Lv;                                                      % Latent heat of vaporization [J / mol]
    boltz = CONSTANTS.boltz;                                                % [W/m^2/K^4] Stefan-Boltzmann constant 
    cp_mol = CONSTANTS.cp_mol;                                              % [J/mol/K] specific heat of air at constant pressure  
    R = CONSTANTS.R;                                                        % [J mol^-1 K^-1] Gas Constant
    rho_dry_air = CONSTANTS.rho_dry_air;                                    % [kg/m3] Density Dry Air
    mmH2OtoMPa = CONSTANTS.mmH2OtoMPa;                                      % Conversion mm to MPa
    
    %SWITCHES
    useG_on=SWITCHES.useG_on;                                               % Use of G measured in the soil energy balance instead of compute it
    useTs_on=SWITCHES.useTs_on;                                             % Use of Ts (soil temperature) [soil energy balance]   
    
    %FORCING
    if (useTs_on && ~isnan(FORCING.Ts(1))) 
        Ts1(1)=FORCING.Ts(1);                                               % Allocate Soil Temperature First Layer
        Ts1(2)=FORCING.zTs(1);                                              % Allocate Node Depth First Layer
    else
        Ts1 = VARIABLES.SOIL.Ts(1);                                         % Allocate Soil Temperature First Layer 
    end
    if useG_on
       G=FORCING.G;                                                         % Allocate Ground Heat Flux        
    else
        G=nan;
    end
    
    
%*************************************************************************
%*************************************************************************


psis1_MPa = psis1 * mmH2OtoMPa;

% Solve SEB_Remainder with fzero
Ts = fzero(@(Ts) SEB_Remainder(Ts, Rabs, Ta1, ea1, pa1, U1, z1, dzs1,...
psis1_MPa, vonk, z0, TC1, epss,Ts1,G,SWITCHES), Ta1);

% Soil Surface Energy Balance
    rhoa = rho_dry_air  * 1000 / 28.97; % [mol/m^3]
    cp = 29.3;                          % specific heat of air at constant pressure [J / mol / K]
    Vw = 18;                            % [g/mol]
    
    
%Gs = Rabs - LEs - Hs;
    

    % Heat and vapor exchange coefficient (Eqn. 9)
    D = U1*vonk^2 / (log(z1/z0))^2;
    
    if useTs_on
       Gs = TC1 * (Ts-Ts1) / dzs1 / 2 ;              
    else
        Gs = TC1 * (Ts-Ts1) / dzs1 / 2;        
    end    
    
    if useG_on
        Gs=G;
    end
    
    Hs = cp*rhoa*D*(Ts-Ta1);    
    esatTs = 0.611*exp(17.502*Ts/(Ts+240.97));
    RH = exp(psis1_MPa*Vw/R/(Ts+273.15));
    LEs = Lv*rhoa*D*(0.622/pa1)*(esatTs * RH - ea1);
        
    
    LWups = epss*boltz*(Ts+273.15)^4;
    
    
    remain = Rabs - Hs - LEs - Gs - LWups;  
 %    remain = Rnrad_soil - Hs - LEs - Gs;  
    
    VARIABLES.SOIL.LEl = 0;                                                 % Latent Heat from Snow-Litter Pack [W/m2]
    VARIABLES.SOIL.LEs = LEs;                                               % Latent Heat from The Soil [W/m2]    
    VARIABLES.SOIL.Hl = 0;                                                  % Sensible heat from Snow-Litter pack [W/m2]  
    VARIABLES.SOIL.Hs = Hs;                                                 % Sensible heat from Snow-Litter pack [W/m2]  
    
    VARIABLES.SOIL.RHl = nan;                                               % Relative Humidity in snow-litter pack.
    VARIABLES.SOIL.Tli = Ts;                                                % Temperature Snow-Litter pack [C]  
    VARIABLES.SOIL.Tsl = Ts;                                                % Tempeature in the Soil-Litter Interface [C]
    VARIABLES.SOIL.Gsl = 0;                                                 % Ground Heat Flux into the Soil - Snow-Litter Boundary from The Snow-Litter Pack [W/m2]
    VARIABLES.SOIL.Gs = Gs;                                                 % Ground Heat Flux into the Soil from the Soil- Snow-Litter Pack  Boundary [W/m2]
    VARIABLES.SOIL.psili = nan;                                             % Soil Water Potential in the snow-Litter Pack [converted from m to mm]  
    VARIABLES.SOIL.psil_MPa = nan;                                          % Soil Water Potential in the snow-Litter Pack [MPa]
    VARIABLES.SOIL.remain = remain;                                         % Total Error in the Soil Energy Balance [W/m2]
    
    VARIABLES.SOIL.deltawice = 0;                                           % Delta of ice needed for change of phase in the snow-litter pack [W/m2]
    VARIABLES.SOIL.dH = 0;                                                  % Delta of Energy [W/m2]
    VARIABLES.SOIL.dS = 0;                                                  % Delta of Energy [W/m2]
    
    VARIABLES.SOIL.case1sl = nan;                                           % Case 
 