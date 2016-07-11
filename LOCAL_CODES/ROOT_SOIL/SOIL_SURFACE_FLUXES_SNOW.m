function [Hl, LE, Gs, Gsl, RHl, Tsurf, remain, dH, dS, VARIABLES] =...
    SOIL_SURFACE_FLUXES_SNOW(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES, FORCING)


%=========================================================================
% Solve surface energy balance if Snow is Present (No-Litter)
%
% Written by Juan Quijano, UIUC, 2013
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
%       Hl              % [W/m2] Sensible heat from Snow pack (or soil)   
%       LE              % [W/m2] Latent heat from Snow pack (or soil)
%       Gs              % [W/m2] Ground heat flux   
%       Gsl             % [W/m2] Ground Heat Flux into the Soil-Snow Boundary    
%       RHl             % [] Relative Humidity of Snow Pack  
%       Tsurf           % [C] Temperature in the surface (Soil-Snow Interface)
%       Remain          % [W/m2] Error in soil energy balance
%       dH              % Delta of Energy Melting/Fusion
%       dS              % Delta of Energy in Time
%       VARIABLES       % VARIABLES structure
% 
%========================================================================


%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % CONTANTS
    dtime = CONSTANTS.dtime;                                                %[s] Time Step 
    vonk = CONSTANTS.vonk;                                                  %[] von Karman constant
    rho_dry_air = CONSTANTS.rho_dry_air;                                    %[kg/m3] Density Dry Air
    mmH2OtoMPa = CONSTANTS.mmH2OtoMPa;                                      % Conversion mm to MPa       
    cp = CONSTANTS.cp_mol;                                                  %[J / mol / K] specific heat of air at constant pressure 
    Lv = CONSTANTS.Lv;                                                      %[J / mol] latent heat of vaporization 
    Lf = CONSTANTS.Lf;                                                      %[J / mol] latent heat of fusion 
    Lf_kg = CONSTANTS.Lf_kg;                                                %[J / kg] latent heat of fusion   
    boltz = CONSTANTS.boltz;                                                %[W/m^2/K^4] Stefan-Boltzmann constant 
    

    % VARIABLES structure
    Ta1 = VARIABLES.CANOPY.TAz(1);                                          %[C] Temperature in the Atmosphere first layer
    ea1 = VARIABLES.CANOPY.EAz(1);                                          %[kPa] Vapor pressure at surface
    pa1 = VARIABLES.CANOPY.PAz(1);                                          %[kPa] Air pressure at surface
    U1 = VARIABLES.CANOPY.Uz(1);                                            %[m/s] wind speed at surface 
    Rabs = VARIABLES.SOIL.Totabs_soil;                                      %[W/m2] Radiation Energy Absorbed Soil 
    psis1 = VARIABLES.SOIL.smp(1);                                          %[mm] Soil water potential first layer
    TC1 = VARIABLES.SOIL.TKsoil(1);                                         %[W / m / K] Thermal Conductivity firs soil layer    
    Ts1 = VARIABLES.SOIL.Ts(1);                                             %[C] Temperature in the Soil fist layer  
    wliqsl = VARIABLES.SOIL.wliqsl;                                         %[kg / m^2] Liquid Water Density per unit area
    wicesl = VARIABLES.SOIL.wicesl;                                         %[kg / m^2] ice Water Density per unit area
    rhosn = VARIABLES.SOIL.rhosn;                                           %[kg / m^3] Snow Density  
    zsn = VARIABLES.SOIL.zsn;                                               %[mm] Snow Depth 
    Tsl =  VARIABLES.SOIL.Tsl;                                              %[C] Top Soil Temperature 
    Tlprev = VARIABLES.SOIL.Tlprev;                                         %[C] Previous time step snow temperature 
    
    
    % VERTSTRUCT structure
    z1 = VERTSTRUC.znc(1);                                                  % [m] Height of bottom atmosphere node
    dzs1 = VERTSTRUC.dzs(1);                                                % [m] Thickness of top soil layer
    
    % PARAMETER structure        
    epss = PARAMS.Rad.epss;                                                 % epss Soil emissivity
    z0 = PARAMS.Soil.z0;                                                    % [m] surface roughness length 
    TK_ice = PARAMS.Soil.TK_ice;                                            % [W / m / K] Thermal Conductivity Ice
    TK_air = PARAMS.Soil.TK_air;                                            % [W / m / K] Thermal Conductivity Air
    HC_liq = PARAMS.Soil.HC_liq;                                            % [J / kg / K] Heat Capacity Liq    
    HC_ice = PARAMS.Soil.HC_ice;                                            % [J / kg / K] Heat Capacity Ice  
    thetamin = PARAMS.Soil.thetamin;                                        % [] Value of soil moisture at which evaporation becomes negligible
    kksnow = PARAMS.Soil.kksnow;                                            % [1/cm] Radiation attenuation of snow

    psis1_MPa = psis1 * mmH2OtoMPa;   
%*************************************************************************
%  Calculate the thermal properties in the snow-litter pack and the 2 cases
%  i. Calculate Thermal Conductivity of snow pack
        TK_sl = TK_air+(7.75*10^(-5)*rhosn+1.105*10^(-6)*rhosn^2)*(TK_ice-TK_air);
% ii. Calculate the volumetric specific heat     
        VHC_SL =  wicesl/zsn*HC_ice + wliqsl/zsn*HC_liq;
% iii. Calculate the variable of both cases.
        zsl_m = zsn/1000;
        kksl = kksnow;

%*************************************************************************    
% calculate the fraction of solar radiation that is absorbed at the
% different depths
Af1 = 1 - exp(-kksl*zsl_m*100);                 % Fraction in the snow layer
Af2 = 1 -  Af1;                                 % Fraction in the soil - snow interface
Af1 = 1;
Af2 = 0;
%*************************************************************************
%  Calculate the snow-litter pack energy balance by an iterative approach
            
dif = 1000;
count = 0;
dH = 0;
Tli = 0;
while abs(dif) > 0.005
         count = count + 1;
         if count > 10
             stop
         end
       dH = 0;                                           
        %i. Compute the energy balance in the surface of the snow-litter pack
       [remain1, Hl, LEl, Gsl, ra, LWups] = SEB_Remainder_snow(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                   Rabs, dH, Ta1, ea1, pa1, U1, z1,...
                                   vonk, z0, epss, zsl_m, TK_sl,...
                                   rho_dry_air, cp, Lv, Lf, boltz, thetamin);                                  
                                   

                                
                               
        Hsl = Rabs*Af1 - Gsl - LEl - Hl - LWups - dH - (VHC_SL*(zsl_m)/dtime)*(Tli-Tlprev); % [W/m2]
        Hmsl = Hsl*dtime/Lf_kg;                                                                % [kg/m2] of ice/water
        if (wicesl > 0 && Hmsl > 0)
            wiceslnew = max(wicesl - Hmsl , 0);
            indexss = Hmsl > wicesl;                                                           % Index to check if snow smelt completely
        elseif (wliqsl > 0 && Hmsl < 0)   
            wiceslnew = min(wicesl - Hmsl, wicesl + wliqsl);
            indexss = logical(0);            
        else
            wiceslnew = wicesl;
            indexss = logical(0);
        end
        deltawice = wiceslnew - wicesl;
        Hstar = Hsl + Lf_kg*deltawice/dtime;
        dH = -Lf_kg*deltawice/dtime;
        Tli = fzero(@(Tli) SEB_Remainder_snow(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                   Rabs, dH, Ta1, ea1, pa1, U1, z1,...
                                   vonk, z0, epss, zsl_m, TK_sl,...
                                   rho_dry_air, cp, Lv, Lf, boltz, thetamin), Ta1);

       [remain1, Hl, LEl, Gsl, ra, LWups] = SEB_Remainder_snow(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                   Rabs, dH, Ta1, ea1, pa1, U1, z1,...
                                   vonk, z0, epss, zsl_m, TK_sl,...
                                   rho_dry_air, cp, Lv, Lf, boltz, thetamin);                                  
                               
        %Tli = (dtime/(VHC_SL*(zsl_m)))*Hstar;
%        if abs(Tli>1e-4)
%           Tli = 0;
%        end
        %ii. solve the energy balance in the bottom of snow litter pack.
                                   
        [remain2, LEs, Gsl, Gs, Tsl] = SEB_Remainder_soil3(Tli, Ts1,...
                 zsl_m, TK_sl, dzs1, TC1, Af2*Rabs);                 
                               
         if count > 1
            dif1 = ((Tli+273) - (Tliprev+273))/(Tliprev+273);                                   
            dif2 = abs(dH - dHprev); 
            dif = max(dif1,dif2);
            if indexss == 1
                dH = wicesl*Lf_kg/dtime;
                wiceslnew = 0;
                deltawice = wiceslnew - wicesl;                
                break;
            end                
         end 
         dHprev= dH;         
         Tliprev = Tli;         
end        
dS = (VHC_SL*(zsl_m)/dtime)*(Tli-Tlprev);
remain = (remain1 + remain2);
LE = LEs + LEl;
Tsurf = Tli;
Hs = 0;
RHl = nan;
%**************************************************************************
% If all Snow melted, then solve energy balance only in the soil surface
if indexss     
    useG_on = SWITCHES.useG_on;          % Use of G measured in the soil energy balance instead of compute it
    if useG_on
       G=FORCING.G;        
    else
        G=nan;
    end
    % Solve SEB_Remainder with fzero
    
    Rabs_mod = Rabs - dH - dS;
    [Ts] = fzero(@(Ts) SEB_Remainder(Ts, Rabs_mod, Ta1, ea1, pa1, U1, z1, dzs1,...
    psis1_MPa, vonk, z0, TC1, epss,Ts1,G,SWITCHES), Ta1);
    
    [remain, Hs, LEs, Gs, Ts, LWups] = SEB_Remainder(Ts, Rabs_mod, Ta1, ea1, pa1, U1, z1, dzs1,...
    psis1_MPa, vonk, z0, TC1, epss,Ts1,G,SWITCHES);
     LEl = 0;
     LE = LEs + LEl;     
     Hl = 0;
     RHl = nan;
     Tli = Ts;
     Tsl = Ts;
     Gsl = 0;
     psil = nan;
     Hstar = dH;
     Hsl = dH;
     Hmsl = wicesl;
     Tsurf = Ts;
end    
    % Save Variables
       
    VARIABLES.SOIL.LEl = LEl;              % Latent Heat from Snow-Litter Pack [W/m2]
    VARIABLES.SOIL.LEs = LEs;              % Latent Heat from The Soil [W/m2]    
    VARIABLES.SOIL.Hl = Hl;                % Sensible heat from Snow-Litter pack [W/m2]  
    VARIABLES.SOIL.Hs = Hs;                % Sensible heat from Snow-Litter pack [W/m2]  
    
    VARIABLES.SOIL.RHl = nan;              % Relative Humidity in snow-litter pack.
    VARIABLES.SOIL.Tli = Tli;              % Temperature Snow-Litter pack [C]  
    VARIABLES.SOIL.Tsl = Tsl;              % Tempeature in the Soil-Litter Interface
    VARIABLES.SOIL.Gsl = Gsl;              % Ground Heat Flux into the Soil - Snow-Litter Boundary from The Snow-Litter Pack [W/m2]
    VARIABLES.SOIL.Gs = Gs;                % Ground Heat Flux into the Soil from the Soil- Snow-Litter Pack  Boundary [W/m2]
    VARIABLES.SOIL.psili = nan;            % Soil Water Potential in the snow-Litter Pack [converted from m to mm]  
    VARIABLES.SOIL.psil_MPa = nan;         % Soil Water Potential in the snow-Litter Pack [MPa]
    VARIABLES.SOIL.remain = remain;        % Total Error in the Energy Balance
    
    VARIABLES.SOIL.deltawice = deltawice;  % Delta of ice needed to change phase in the snow-litter pack
    VARIABLES.SOIL.dH = dH;                % Delta of Energy Melting/Fusion 
    VARIABLES.SOIL.dS = dS;                % Delta of Energy in Time 
    
    VARIABLES.SOIL.case1sl = nan;          % Case 
