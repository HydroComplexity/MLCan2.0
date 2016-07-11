function [Hl, LE, Gs, Gsl, RHl, Tli, remain, dH, dS, VARIABLES] =...
    SOIL_SURFACE_FLUXES_LITTER(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES, FORCING)


%=========================================================================
% Solve surface energy balance including a Litter Layer
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
%       Hl              % [W/m2] Sensible heat from Snow-Litter pack   
%       LE              % [W/m2] Latent heat from Snow-Litter pack   
%       Gs              % [W/m2] Ground heat flux   
%       Gsl             % [W/m2] Ground Heat Flux into the Soil - Snow-Litter Boundary    
%       RHl             % [] Relative Humidity of Snow-Litter Pack  
%       Tli             % [C] Temperature in the surface (Soil-Litter Interface)
%       Remain          % [W/m2] Error in soil energy balance
%       dH              % Delta of Energy Melting/Fusion
%       dS              % Delta of Energy in Time
%       VARIABLES       % VARIABLES structure
% 
%========================================================================


%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % CONSTANTS
    dtime = CONSTANTS.dtime;                                                % [s] Time Step
    vonk = CONSTANTS.vonk;                                                  % von Karman constant
    rho_dry_air = CONSTANTS.rho_dry_air;                                    % %[kg/m3] Density Dry Air 
    mmH2OtoMPa = CONSTANTS.mmH2OtoMPa;                                      % Conversion from mm to MPa
    cp = CONSTANTS.cp_mol;                                                  % [J / mol / K] specific heat of air at constant pressure 
    Lv = CONSTANTS.Lv;                                                      % [J / mol] Latent heat of vaporization 
    Lf = CONSTANTS.Lf;                                                      % [J / mol] Latent heat of fusion 
    Lf_kg = CONSTANTS.Lf_kg;                                                % [J / kg] Latent heat of fusion    
    boltz = CONSTANTS.boltz;                                                % [W/m^2/K^4] Stefan-Boltzmann constant
    Vw = CONSTANTS.Vw;                                                      % Mols to kg conversion 
    R = CONSTANTS.R;                                                        % [J mol^-1 K^-1] Gas Constant
    
    % VARIABLES
    Ta1 = VARIABLES.CANOPY.TAz(1);                                          % [C] Atmospheric Temperature
    if abs(Ta1) > 70
        Ta1 = FORCING.Ta;
    end
    ea1 = VARIABLES.CANOPY.EAz(1);                                          % [kPa] Vapor pressure at surface
    pa1 = VARIABLES.CANOPY.PAz(1);                                          % [kPa] Air pressure at surface
    U1 = VARIABLES.CANOPY.Uz(1);                                            % [m/s] wind speed at surface
    Rabs = VARIABLES.SOIL.Totabs_soil;                                      % [W/m2] Radiation Energy Absorbed Soil
    volliq1 = VARIABLES.SOIL.volliq(1);                                     % [] Volumetric Water Content Top Soil
    psis1 = VARIABLES.SOIL.smp(1);                                          % [mm] Soil water potential of top layer
    TC1 = VARIABLES.SOIL.TKsoil(1);                                         % [W / m / K] Thermal Conductivity firs soil layer
    Ts1 = VARIABLES.SOIL.Ts(1);                                             % [C] Temperature in the Soil fist layer 
    dzlit = VARIABLES.SOIL.dzlit_m;                                         % [m] Thickness of litter layer
    dzlit_mm = VARIABLES.SOIL.dzlit_m*1000;                                 % [mm] Thickness of litter layer    
    wliqsl = VARIABLES.SOIL.wliqsl;                                         % [kg / m^2] Liquid Water Density per unit area    
    wicesl = VARIABLES.SOIL.wicesl;                                         % [kg / m^2] Ice Water Density per unit area
    zliqsl = VARIABLES.SOIL.zliqsl;                                         % [mm] Liquid water depth in snow litter pack
    zicesl = VARIABLES.SOIL.zicesl;                                         % [mm] Ice water depth in snow litter pack
    volliqli = VARIABLES.SOIL.volliqli;                                     % [] Volumetric Water Content Litter    
    voliceli = VARIABLES.SOIL.voliceli;                                     % [] Volumetric Ice Content Litter
    rhosn = VARIABLES.SOIL.rhosn;                                           % [kg / m^3] Snow pack density
    zsn = VARIABLES.SOIL.zsn;                                               % [C] Previous Soil Temperature
    Tlprev = VARIABLES.SOIL.Tlprev;
    
    % SWITCHES
    littertype = SWITCHES.littertype;                                       % Type of litter solution 
    
    % VERTSTRUCT structure
    Tsl = VARIABLES.SOIL.Tsl;                                               % [C] Top Soil Temperature
    
    % INITIALIZE structure    
    z1 = VERTSTRUC.znc(1);                                                  % [m] height of bottom atmosphere node
    dzs1 = VERTSTRUC.dzs(1);                                                % [m] Thickness first soil layer 
    porsl1 = VERTSTRUC.porsl(1);                                            % [] Porosity Top Soil Layer
    
    % PARAMETER structure        
    ldif = PARAMS.Soil.ldif;                                                % [m/s] Vapor litter diffusivity  
    epss = PARAMS.Rad.epss;                                                 % Soil emissivity
    z0 = PARAMS.Soil.z0;                                                    % [m] surface roughness length        
    c1 = PARAMS.Soil.c1;                                                    % Parameter to compute litter resistance
    thetatr = PARAMS.Soil.thetatr;                                          % Parameter to compute litter resistance
    psill = PARAMS.Soil.psill;                                              % Parameter to compute litter water potential
    bl = PARAMS.Soil.bl;                                                    % Parameter to compute litter water potential        
    bdl = PARAMS.Soil.bdl;                                                  % Parameter to compute litter water potential
    rhowater = PARAMS.Soil.rhowater;                                        % [kg/m3] Density of liquid water 
    TK_litter = PARAMS.Soil.TK_litter;                                      % [W/m/k] Litter Thermal Conductivity                                 
    VHC_litter = PARAMS.Soil.VHC_litter;                                    % [J/m3/K] Volumetric Heat Capacity of Litter   
    TK_liq = PARAMS.Soil.TK_liq;                                            % [W / m / K] Thermal Conductivity liquid water
    TK_ice = PARAMS.Soil.TK_ice;                                            % [W / m / K] Thermal Conductivity ice
    TK_air = PARAMS.Soil.TK_air;                                            % [W / m / K] Thermal Conductivity air
    HC_air = PARAMS.Soil.HC_air;                                            % [J / kg / K] Heat Capacity Air
    HC_liq = PARAMS.Soil.HC_liq;                                            % [J / kg / K] Heat Capacity Liquid water
    HC_ice = PARAMS.Soil.HC_ice;                                            % [J / kg / K] Heat Capacity ice
    thetamin = PARAMS.Soil.thetamin;                                        % [] Value of volumetric litter content at which evaporation becomes negligible
    thetals = PARAMS.Soil.thetals;                                          % [] litter porosity
    Tf = PARAMS.Soil.Tf;                                                    % [K] Freexing Temperatures of fresh water
    slfactor = PARAMS.Soil.slfactor;                                        % Percentage above which is considered as snow for energy balance
    kklitter = PARAMS.Soil.kklitter;                                        % [1/cm] Radiation attenuation for litter   
    kksnow = PARAMS.Soil.kksnow;                                            % [1/cm] Radiation attenuation for snow

    dzs1 = dzs1(1);        
    psis1_MPa = psis1 * mmH2OtoMPa;   
%*************************************************************************
%  Calculate the thermal properties in the snow-litter pack and the 2 cases
%  i. Calculate Thermal Conductivity of snow-litter pack
    if zsn <= dzlit_mm
        % Calculate Sr
        Sr = (volliqli + voliceli)/thetals;   % litter effective porosity
        % Calculate Ke
        if Tlprev + Tf < Tf
            Ke = Sr;
        else
            Ke = max(log(Sr)+1, 0);
        end
        % compute TKsat
        if Tlprev + Tf < Tf
            TKsat = TK_litter^(1-thetals)*TK_liq^(thetals)*TK_ice^(thetals-volliqli);
        else
            TKsat = TK_litter^(1-thetals)*TK_liq^(thetals);
        end    
        TK_sl = Ke*TKsat+(1-Ke)*TK_litter;
    else
        TK_sl = TK_air+(7.75*10^(-5)*rhosn+1.105*10^(-6)*rhosn^2)*(TK_ice-TK_air);
    end
% ii. Calculate the volumetric specific heat     
    if zsn <= dzlit_mm
        VHC_SL = VHC_litter*(1-thetals) + wicesl/dzlit*HC_ice + wliqsl/dzlit*HC_liq;        
    else
        %VHC_SL = VHC_litter*(1-thetals)*dzlit/zsn + wicesl/zsn*HC_ice + wliqsl/zsn*HC_liq;
        VHC_SL =  wicesl/zsn*HC_ice + wliqsl/zsn*HC_liq;
    end        
% iii. Calculate the variable of both cases.
    if zicesl < slfactor*dzlit_mm*thetals
        case1sl = false;                               % Assumed as No ice   
    else
        case1sl = true;                                % Ice
    end
    
% Define the size of zsl_m     
    if zsn < dzlit_mm*thetals                       % Lower than the litter 
        zsl_m = dzlit;
        kksl = kklitter;
    else
        zsl_m = zsn/1000 + dzlit*(1-thetals);       % Higher than the litter
        kksl = kksnow;
    end

%*************************************************************************    
% calculate the fraction of solar radiation that is absorbed at the
% different depths
Af1 = 1 - exp(-kksl*zsl_m*100);                 % Fraction in the snow layer
Af2 = 1 -  Af1;                                 % Fraction in the soil - snow interface
Af1 = 1;
Af2 = 0;
%*************************************************************************
%  Calculate the snow-litter pack energy balance by an iterative approach
            
%criteria = 1000;
dif = 1000;
count = 0;
Tli = 0 ;                   % Start assuming 0 to get maximum energy to freeze or melt
while abs(dif) > 0.005
         count = count + 1;
         if count > 10
             stop =34;
         end
        
        %i. Compute the energy balance in the surface of the snow-litter pack
       dH = 0;                                           
        [remain1, Hl, LEl, Gsl, RHl, ra, psil_MPa, psil, LWups] = SEB_Remainder_litter(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                       Af1*Rabs, dH, Ta1, ea1, pa1, U1, z1,...
                                       vonk, z0, epss, c1, thetatr, psill, bl,...
                                       bdl, rhowater, zsl_m, TK_sl, volliqli,...
                                       rho_dry_air, cp, Lv, Lf, boltz, Vw, R, mmH2OtoMPa, thetamin,...
                                       case1sl);                                   
                                   

        Hsl = Rabs*Af1 - Gsl - LEl - Hl - LWups - dH - (VHC_SL*(zsl_m)/dtime)*(Tli-Tlprev); % [W/m2]
        Hmsl = Hsl*dtime/Lf_kg;                                                                % [kg/m2] of ice/water
        if (wicesl > 0 && Hmsl > 0)
            wiceslnew = max(wicesl - Hmsl , 0);
        elseif (wliqsl > 0 && Hmsl < 0)   
            wiceslnew = min(wicesl - Hmsl, wicesl + wliqsl);
        else
            wiceslnew = wicesl;
        end
        deltawice = wiceslnew - wicesl;
        Hstar = Hsl + Lf_kg*deltawice/dtime;
        dH = -Lf_kg*deltawice/dtime;
        Tli = fzero(@(Tli) SEB_Remainder_litter(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                       Af1*Rabs, dH, Ta1, ea1, pa1, U1, z1,...
                                       vonk, z0, epss, c1, thetatr, psill, bl,...
                                       bdl, rhowater, zsl_m, TK_sl, volliqli,...
                                       rho_dry_air, cp, Lv, Lf, boltz, Vw, R, mmH2OtoMPa, thetamin,...
                                       case1sl), Ta1);
        
        [remain1, Hl, LEl, Gsl, RHl, ra, psil_MPa, psil, LWups] = SEB_Remainder_litter(Tli, Tsl, Tlprev, VHC_SL, dtime,...
                                       Af1*Rabs, dH, Ta1, ea1, pa1, U1, z1,...
                                       vonk, z0, epss, c1, thetatr, psill, bl,...
                                       bdl, rhowater, zsl_m, TK_sl, volliqli,...
                                       rho_dry_air, cp, Lv, Lf, boltz, Vw, R, mmH2OtoMPa, thetamin,...
                                       case1sl);                                   
        
        
        %Tli = (dtime/(VHC_SL*(zsl_m)))*Hstar;
        if abs(Tli)<1e-4
           Tli = 0;
        end
        %ii. solve the energy balance in the bottom of snow litter pack.
        %phig1=(1-volliq1)/(porsl1);
        %Dsoil1 = phig1*0.66*sdif;    
                                   
        if littertype == 1
            
             if case1sl == 1
                       [remain2, LEs, Gsl, Gs, Tsl] = SEB_Remainder_soil3(Tli, Ts1,...
                                   zsl_m, TK_sl, dzs1, TC1, Af2*Rabs);   
                       if abs(Gs) > 4*abs(Rabs)
                               Gscorr = 0.05*abs(Rabs)*sign(Gs);
                               [remain2, LEs, Gsl, Gs, Tsl] = SEB_Remainder_soil3_correction(Tli, Ts1,...
                                           zsl_m, TK_sl, dzs1, TC1, Af2*Rabs, Gscorr);   
                       end   
                               
                               
             elseif case1sl == 0  
                          
                 Tsl = fzero(@(Tsl) SEB_Remainder_soil1(Tli, Tsl, Ts1, ea1, pa1, U1, z1,...
                                       vonk, z0, zsl_m, TK_sl,...
                                       rho_dry_air, psis1_MPa, dzs1, TC1, Af2*Rabs), Ta1);

                [remain2, LEs, Gsl, Gs] = SEB_Remainder_soil1(Tli, Tsl, Ts1, ea1, pa1, U1, z1,...
                                       vonk, z0, zsl_m, TK_sl,...
                                       rho_dry_air, psis1_MPa, dzs1, TC1, Af2*Rabs);                                   
             end                                   
                                    
        elseif littertype == 2
            
            if case1sl == 1
                       [remain2, LEs, Gsl, Gs, Tsl] = SEB_Remainder_soil3(Tli, Ts1,...
                                   zsl_m, TK_sl, dzs1, TC1, Af2*Rabs);                 
            else                
            
                        Tsl = fzero(@(Tsl) SEB_Remainder_soil2(Tl, Tsl, Ts1, porsl1,...
                                   pa1, c1, thetatr, psill, bl,...
                                   bdl, rhowater, zsl_m, TK_litter, ldif, volliqli,...
                                   rho_dry_air, mmH2OtoMPa, psis1_MPa, dzs1, TC1, Af2*Rabs), Ta1);

                        [remain2, LEs, Gsl, Gs] = SEB_Remainder_soil2(Tl, Tsl, Ts1, porsl1,...
                                   pa1, c1, thetatr, psill, bl,...
                                   bdl, rhowater, zsl_m, TK_litter, ldif, volliqli,...
                                   rho_dry_air, mmH2OtoMPa, psis1_MPa, dzs1, TC1, Af2*Rabs);                               
            end            
        end   
                               
         if count > 1
            dif1 = ((Tli+273) - (Tliprev+273))/(Tliprev+273); 
            dif2 = abs(dH - dHprev);
            dif = max(dif1,dif2);
%            criteria = abs(dif-difprev)/dif; 
%            difprev = dif;
           %dif = (Tli - Tliprev)/Tliprev;                                   
         end 
         Tliprev = Tli; 
         dHprev= dH;         
end        
    dS = (VHC_SL*(zsl_m)/dtime)*(Tli-Tlprev);
    remain = (remain1 + remain2);
    LE = LEs + LEl;    
    
    % Define the structure of the litter
        
    VARIABLES.SOIL.LEl = LEl;                                              % [W/m2] Latent Heat from Snow-Litter Pack 
    VARIABLES.SOIL.LEs = LEs;                                              % [W/m2] Latent Heat from The Soil     
    VARIABLES.SOIL.Hl = Hl;                                                % [W/m2] Sensible heat from Snow-Litter pack   
    VARIABLES.SOIL.Hs = nan;                                               % [W/m2] Sensible heat from Soil   
    
    VARIABLES.SOIL.RHl = RHl;                                              % Relative Humidity in snow-litter pack.
    VARIABLES.SOIL.Tli = Tli;                                              % [C] Temperature Snow-Litter pack   
    VARIABLES.SOIL.Tsl = Tsl;                                              % [C] Tempeature in the Soil-Litter Interface
    VARIABLES.SOIL.Gsl = Gsl;                                              % [W/m2] Ground Heat Flux into the Soil - Snow-Litter Boundary from The Snow-Litter Pack 
    VARIABLES.SOIL.Gs = Gs;                                                % [W/m2] Ground Heat Flux into the Soil from the Soil- Snow-Litter Pack  Boundary 
    VARIABLES.SOIL.psili = psil*1000;                                      % [mm] Soil Water Potential in the snow-Litter Pack [converted from m to mm]  
    VARIABLES.SOIL.psil_MPa = psil_MPa;                                    % [MPa] Soil Water Potential in the snow-Litter Pack 
    VARIABLES.SOIL.remain = remain;                                        % [W/m2] Total Error in the Energy Balance
    
    VARIABLES.SOIL.deltawice = deltawice;                                  % [] Delta of ice needed for change of phase in the snow-litter pack
    VARIABLES.SOIL.dH = dH;                                                % [W/m2] Delta of Energy Melting/Fusion
    VARIABLES.SOIL.dS = dS;                                                % [W/m2] Change of Internal Energy in Snow-Litter Pack in Time Step

    VARIABLES.SOIL.case1sl = case1sl;                                      % Case for snow solution
    
