function [remain, Hl, LEl, G1, ra, LWups] = SEB_Remainder_snow(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                   Rabs, dH, Ta1, ea1, pa, U1, z1,...
                                   vonk, z0, epss, zsl, TK_sl,...
                                   rho_dry_air, cp, Lv, Lf, boltz, thetamin)

                              
                              
%=========================================================================
% Solve surface energy in a snow pack.
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       Tli             % [C] Temperature in the litter layer 
%       Tsl             % [C] Temperature in the Litter Snow - Soil interface 
%       Tlprev          % [C] Temperature in the litter-snow pack in the previous time step 
%       VHC_SL          % [J / m^3 / K] volumetric specific heat 
%       dtime           % [s] Time Step
%       Rabs            % [W/m2] Energy Absorption Soil
%       dH              % [W/m2] Energy required/release for melting/freezing
%       Ta1             % [C] Temperature in the Atmosphere first layer 
%       ea1             % [kPa] Vapor pressure at surface 
%       pa              % [kPa] Air pressure at surface [kPa]
%       U1              % [m/s] wind speed at surface 
%       z1              % [m] height of bottom atmosphere node 
%       vonk            % [] von Karman constant 
%       z0              % [m] surface roughness length 
%       epss            % Soil emissivity
%       zsl             % [m] Depth snow pack 
%       TK_sl           % [W / m / K] Thermal Conductivity of snow-litter pack 
%       rho_dry_air     % [kg/m3] Density Dry Air 
%       cp              % specific heat of air at constant pressure 
%       Lv              % Latent heat of vaporization [J / mol] 
%       Lf              % Latent heat of vaporization [J / mol] 
%       boltz           % [W/m^2/K^4] Stefan-Boltzmann constant  
%       thetamin        % [] Value of soil moisture at which evaporation becomes negligible     
%------------------------- Output Variables ------------------------------
%       remain          % [W/m2] Energy Balance Error   
%       Hl              % [W/m2] Sensible heat from Snow Litter pack    
%       LEl             % [W/m2] Latent heat from Snow Litter pack    
%       Gl              % [W/m2] Ground Heat Flux into the Snow Litter pack
%       ra              % [s / m] Resistivity for water vapor  Ogee and Brunet [2002]
%       LWups           % [W / m2] Emmited Longwave Radiation
%========================================================================  

% Soil Surface Energy Balance
    rhoa = rho_dry_air  * 1000 / 28.97; % from [kg/m3] to [mol/m^3]    
    Ls = Lv + Lf;

    % COMPUTATION OF LATENT HEAT FROM LITTER
    % compute aerodynamic resistance 
    ra = ((log(z1/z0))^2)/(U1*vonk^2); %  D = U1*vonk^2 / (log(z1/z0))^2;
    % Compute saturated pressure at the litter temperature
    esatTl = 0.611*exp(17.502*Tli/(Tli+240.97));
    %compute LE
        LEl = (Ls*rhoa*(0.622/pa)*(esatTl - ea1))/(ra);   
    % COMPUTATION OF SENSIBLE HEAT FROM LITTER
    Hl = cp*rhoa*(Tli-Ta1)/ra;    
    
    % COMPUTATION OF HEAT INTO THE LITTER
    G1 = (TK_sl)*(Tli - Tsl)/((zsl)/2);
    
    % COMPUTATION OF STORAGE TERM
    dS = (VHC_SL*(Tli-Tlprev)/dtime)*zsl;
    
    % COMPUTE REMAIN        
    
    LWups = epss*boltz*(Tli+273.15)^4;
    remain = Rabs - Hl - LEl - G1 - LWups - dH - dS;  
 %    remain = Rnrad_soil - Hs - LEs - Gs;  
 