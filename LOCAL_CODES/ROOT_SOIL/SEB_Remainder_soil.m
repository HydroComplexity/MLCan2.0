function [Remain, LEs, Gsl, Gs] = SEB_Remainder_soil(Tl, Tsl, Ts1, porsl1, Ta1, ea1,...
                                   pa, U1, z1, vonk, z0, c1, thetatr, psill, bl,...
                                   bdl, rhowater, litterthickness, TK_litter, ldif, Dsoil1, volliql,...
                                   rho_dry_air, mmH2OtoMPa, psis1_MPa, dzs1, TC1, ra)
                              
% Solve surface energy balance based on formulation of (Hinzman et al, JGR 1998)
%   --> Equation #s refer to Hinzman
%
%   INPUTS:
%       Ts = soil surface temp [C]
%       Rabs = absorbed radiation by soil [W/m^2]
%       Ta1 = air temperature at surface [C]
%       ea1 = vapor pressure at surface [kPa]
%       pa = air pressure at surface [kPa]
%       U1 = wind speed at surface [m/s]
%       z1 = height of bottom atmosphere node [m]
%       dzs1 = thickness of top soil layer [m]
%       psis1_MPa = soil water potential of top layer [MPa] --> 1MPa = 1J/g
%       vonk = von Karman constant
%       z0 = surface roughness length
%       TC1 = thermal conductivity of top soil layer 

%ldif = ldif * 1000000;  % Change from m2/s to mm2/s
%Dsoil1 = Dsoil1 * 1000000; % Change from m2/s to mm2/s

% Soil Surface Energy Balance
    rhoa = rho_dry_air  * 1000 / 28.97; % from [kg/m3] to [mol/m^3]
    rhow = 1000 ;   % water density [kg /m3]
    cp = 29.3;       % specific heat of air at constant pressure [J / mol / K]
    Lv_g = 2260*1000; % Lv in [J / (mm m2)]  % latent heat of vaporization 
    boltz = 5.6697 * 10^-8;  % [W m^-2 K^-4]
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
%    LEs2 = (Lv*rhoa*(0.622/pa)*(esatTs1 * RHs1 - ea1))/(ra);   
%    Es_m = (1/((1/ldif)+(1/Dsoil1)))*rhow*(0.622/pa)*(esatTs1 * RHs1 - esatTl * RHl)/((litterthickness + dzs1)/(2));   
%    Es_m = (ldif)*rho_dry_air*(0.622/pa)*(esatTs1 * RHs1 - esatTl * RHl)/((litterthickness + dzs1)/(2));   
    Es_m = (ldif)*rho_dry_air*(0.622/pa)*(esatTs1 * RHs1 - esatTsl * RHsl)/((litterthickness)/(2));      
    Es_mm= Es_m*1000;
    LEs = Es_mm*Lv_g;
    
    % COMPUTE THE GROUND HEAT FLUX FROM THE LITTER TO THE SOIL LITTER
    % SURFACE
    
    Gsl = (TK_litter)*(Tl - Tsl)/((litterthickness)/2); 

    Gs = (TC1)*(Tsl - Ts1)/((dzs1)/2); 
    
    Remain = Gsl - LEs - Gs;
    