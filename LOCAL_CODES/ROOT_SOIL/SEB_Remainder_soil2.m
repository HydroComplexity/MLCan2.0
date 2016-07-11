function [remain, LEs, Gsl, Gs] = SEB_Remainder_soil2(Tl, Tsl, Ts1, porsl1,...
                                   pa, c1, thetatr, psill, bl,...
                                   bdl, rhowater, zsl, TK_litter, ldif, volliql,...
                                   rho_dry_air, mmH2OtoMPa, psis1_MPa, dzs1, TC1, Rabs)
                              
%=========================================================================
% Solve surface energy balance with a snow-litter pack CASE2. (No Latent Heat)
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       Tl              % [C] Temperature in the Snow-Litter pack 
%       Tsl             % [C] Temperature in the Litter Snow - Soil interface 
%       Ts1             % [C] Temperature in the Soil fist layer
%       porsl1          % [] Porosity first layer 
%       pa              % [kPa] Air pressure at surface [kPa]
%       cl              % Parameter to compute psi litter
%       thetatr         % Parameter to compute psi litter 
%       psill           % [m] Parameter to compute psi litter
%       b1              % [] Parameter to compute psi litter
%       rhowater        % [] Parameter to compute psi litter
%       zsl             % [m] Snow-Litter pack thickness 
%       TK_litter       % [W / m / K] Thermal Conductivity Snow Pack 
%       ldif            % [m/s] Vapor litter diffusivity 
%       volliql         % [] Litter volumetric water content
%       rho_dry_air     % [kg/m3] Density Dry Air 
%       mmH2OtoMPa      % Conversion factor from mmH2O to MPa
%       psis1_MPa       % [MPa] Soil water potential first layer
%       dzs1            % [m] Thickness first soil layer 
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
    Lv_g = 2260*1000; % Lv in [J / (mm m2)]  % latent heat of vaporization 
    Vw = 18; % [g/mol]
    R = 8.3143; % [J/mol/K]

    % COMPUTATION OF LATENT HEAT FROM LITTER
    % compute litter resistance
    val = c1*(volliql - thetatr);
    rl = max(0, val);
    % compute soil resistance
    rs = 4104*(porsl1-volliql)-805;  
    if rs<0
        rs=1;
    end
    % psi litter
    psil = - psill*((rhowater/bdl)*volliql)^(-bl); % [m] 
    psil_MPa = psil*1000* mmH2OtoMPa;
    % compute relative humidity in the litter layer assuming phase
    % equilibrium
    RHsl = exp(psil_MPa*Vw/R/(Tsl+273.15));
    % compute relative humidity in the soil layer assuming phase
    % equilibrium
    RHs1 = exp(psis1_MPa*Vw/R/(Ts1+273.15));
    % Compute saturated pressure at the litter and soil temperature
    esatTsl = 0.611*exp(17.502*Tsl/(Tsl+240.97));
    esatTs1 = 0.611*exp(17.502*Ts1/(Ts1+240.97));
    %compute LE
    Es_m = (ldif)*rho_dry_air*(0.622/pa)*(esatTs1 * RHs1 - esatTsl * RHsl)/((zsl)/(2));      
    Es_mm= Es_m*1000;
    LEs = Es_mm*Lv_g;
    
    % COMPUTE THE GROUND HEAT FLUX FROM THE LITTER TO THE SOIL LITTER
    % SURFACE
    
    Gsl = (TK_litter)*(Tl - Tsl)/((zsl)/2); 

    Gs = (TC1)*(Tsl - Ts1)/((dzs1)/2); 
    
    remain = Gsl - LEs - Gs + Rabs;
    