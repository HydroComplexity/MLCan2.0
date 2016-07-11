function [remain, LEs, Gsl, Gs] = SEB_Remainder_soil1(Tl, Tsl, Ts1, ea1, pa, U1, z1,...
                                   vonk, z0, zsl, TK_sl,...
                                   rho_dry_air, psis1_MPa, dzs1, TC1, Rabs)

                              
                              
%=========================================================================
% Solve surface energy balance with a snow-litter pack CASE1. (No Latent Heat)
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       Tl              % [C] Temperature in the Snow-Litter pack 
%       Tsl             % [C] Temperature in the Litter Snow - Soil pack 
%       Ts1             % [C] Temperature in the Soil fist layer pack 
%       eal             % [kPa] Vapor pressure at surface 
%       Pal             % [kPa] Air pressure at surface [kPa]
%       Ul              % [m/s] wind speed at surface
%       zl              % [m] Height of bottom atmosphere node [m] 
%       vonk            % [] von Karman constant 
%       z0              % [m] surface roughness length 
%       zsl             % [m] Snow-Litter pack thickness 
%       TK_sl           % [W / m / K] Thermal Conductivity Snow-Litter Pack
%       rho_dry_air     % [kg/m3] Density Dry Air 
%       psis1_MPa       % [MPa] Soil water potential first layer
%       dzs1            % [m] Tickness First Layer Soil  
%       TC1             % [W / m / K] Thermal Conductivity firs soil layer
%       Rabs            % [W/m2] Absorbed Radiation
%------------------------- Output Variables ------------------------------
%       remain          % [W/m2] Energy Balance Error   
%       LEs             % [W/m2] Latent heat from soil (should be zero)   
%       Gs1             % [W/m2] Ground heat flux into the interface   
%       Gs              % [W/m2] Ground Heat Flux into the Soil     
% 
%========================================================================                               


% Soil Surface Energy Balance
    rhoa = rho_dry_air  * 1000 / 28.97; % from [kg/m3] to [mol/m^3]
    Lv = 44000;  % latent heat of vaporization [J / mol]
    Vw = 18; % [g/mol]
    R = 8.3143; % [J/mol/K]

    %   COMPUTE LATENT HEAT FROM SOIL TO ATMOSPHERE
        % COMPUTATION OF LATENT HEAT FROM LITTER
    % compute aerodynamic resistance [s/m]
    ra = ((log(z1/z0))^2)/(U1*vonk^2); %  D = U1*vonk^2 / (log(z1/z0))^2;

    esatTsl = 0.611*exp(17.502*Tsl/(Tsl+240.97));
    RH = exp(psis1_MPa*Vw/R/(Tsl+273.15));

    %compute LE
    LEs = (Lv*rhoa*(0.622/pa)*(esatTsl * RH - ea1))/(ra);   
    
    % COMPUTE THE GROUND HEAT FLUX FROM THE LITTER TO THE SOIL LITTER
    % SURFACE
    
    Gsl = (TK_sl)*(Tl - Tsl)/((zsl)/2); 

    Gs = (TC1)*(Tsl - Ts1)/((dzs1)/2); 
    
    remain = Gsl - LEs - Gs + Rabs;
    