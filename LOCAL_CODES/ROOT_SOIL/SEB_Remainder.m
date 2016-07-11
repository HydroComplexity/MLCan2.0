function [remain, Hs, LEs, Gs, Ts, LWups] = SEB_Remainder(Ts, Rabs, Ta1, ea1, pa, U1, z1, dzs1, ...
                                  psis1_MPa, vonk, z0, TC1, epss , Ts1,G,SWITCHES)

%=========================================================================
% Solve surface energy balance below a snow pack. (No Latent Heat)
%
% Written by Darren Drewry, Modified by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       Ts              % [C] Temperature in the Soil surface 
%       Rabs            % [W/m2] Energy Absorption Soil
%       Ta1             % [C] Temperature in the Atmosphere first layer 
%       ea1             % [kPa] Vapor pressure at surface 
%       pa              % [kPa] Air pressure at surface
%       U1              % [m/s] wind speed at surface 
%       z1              % [m] height of bottom atmosphere node 
%       dzs1            % [m] Thickness of top soil layer 
%       psis1_MPa       % [MPa] Soil water potential of top layer 
%       vonk            % [] von Karman constant 
%       z0              % [m] surface roughness length 
%       TC1             % [W / m / K] Thermal Conductivity firs soil layer
%       epss            % Soil emissivity
%       Ts1             % [C] Soil Temperature first layer
%       G               % [W/m2] Ground Heat Flux 
%------------------------- Output Variables ------------------------------
%       remain          % [W/m2] Energy Balance Error   
%       Ha              % [W/m2] Sensible heat from soil    
%       LEs             % [W/m2] Latent heat from soil    
%       Gs              % [W/m2] Ground Heat Flux into the Soil     
%       Ts              % [C] Soil Surface Temperature 
%       LWups           % [W / m2] Emmited Longwave Radiation
%========================================================================   

    useG_on=SWITCHES.useG_on;          % Use of G measured in the soil energy balance instead of compute it
    useTs_on=SWITCHES.useTs_on;         % Use of Ts (soil temperature) [soil energy balance]   

% Soil Surface Energy Balance
    density_dry_air = 1.2923; % [kg / m^3]
    rhoa = density_dry_air  * 1000 / 28.97; % [mol/m^3]
    cp = 29.3;       % specific heat of air at constant pressure [J / mol / K]
    Lv = 44000;  % latent heat of vaporization [J / mol]
    boltz = 5.6697 * 10^-8;  % [W m^-2 K^-4]
    Vw = 18; % [g/mol]
    R = 8.314; % [J/mol/K]
    
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
    LEs = Lv*rhoa*D*(0.622/pa)*(esatTs * RH - ea1);
    
    
    
    LWups = epss*boltz*(Ts+273.15)^4;
    
    
    remain = Rabs - Hs - LEs - Gs - LWups;  
 %    remain = Rnrad_soil - Hs - LEs - Gs;  
 
