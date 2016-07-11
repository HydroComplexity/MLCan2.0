function [Ph, An, Ci, gsv, Tl, LE, dHcan, TR, Evap_mm, H, psil, fsvg, fsvm, Ch2o_mm, gbv, gbh, ...
          Ph_limit, Jc_C3, Jj_C3, Js_C3, Jc_C4, Jj_C4, Js_C4, ...
          VARIABLES,converged] = ...
    LEAF_SOLUTION(FORCING, VARIABLES, PARAMS, CONSTANTS, VERTSTRUC, SWITCHES, sunlit,cntspecies)
                
%=========================================================================
%   This code solves the leaf dynamics 
%
% Written By: Darren Drewry, Modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING          % FORCING structure
%       VARIABLES        % VARIABLES structure
%       CONSTANTS        % CONSTANTS structure
%       VERTSTRUC        % VERTSTRUC structure
%       SWITCHES         % SWITCHES structure
%       sunlit           % [1] Sunlit [0] Shade
%       cntspecies       % [] Number of species
%------------------------- Output Variables ------------------------------
%       Ph               % [umol CO2/ m^2 leaf / s] Photosynthetic flux from ecosystem 
%       An               % [umol CO2/ m^2 leaf / s] Photosynthetic minus Leaf Respiration flux from ecosystem 
%       gsv              % [mol/m^2/s] Stomatal Conductance 
%       Tl               % [C] Leaf Temperature  
%       LE               % [W / m^2 leaf] Latent Heat From Leaf  
%       dHcan            % [W / m^2 leaf] Heat Storage Leaf
%       TR               % [mm / s / m^2 leaf] Heat Storage Leaf
%       Evap_mm          % [mm / s / m^2 leaf] Evaporation Canopy from 
%       H                % [W / m^2 leaf] Sensible Heat From Leaf
%       psil             % [MPa] Leaf Water Potential 
%       fsvg             % [] g Tuzet factor to Ball Berry Model 
%       fsvm             % [] m Tuzet factor to Ball Berry Model
%       Ch2o_mm          % [mm / s / m^2 leaf] Condensation water in the canopy
%       gbv              % [mol/m^2/s] Boundary Layer conductance to vapor
%       gbh              % [mol/m^2/s] Boundary Layer conductance to heat
%       Ph_limit         % [] Type of Photosynthesis 1. Rubisco-limited, 2light-limited, 3. Sucrose Limited
%       Jc_C3            % [umol/m^2 leaf area/s] Rubisco-Limited Rate for C3
%       Jj_C3            % [umol/m^2 leaf area/s] Light-Limited Rate for C3
%       Js_C3            % [umol/m^2 leaf area/s] Sucrose-Limited Rate for C3
%       Jc_C4            % [umol/m^2 leaf area/s] Rubisco-Limited Rate for C4
%       Jj_C4            % [umol/m^2 leaf area/s] Light-Limited Rate for C4
%       Js_C4            % [umol/m^2 leaf area/s] Sucrose-Limited Rate for C4
%       VARIABLES        % VARIABLES structure
%       converged        % [] Information on the convergence of the leaf subroutine

%========================================================================              
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES
    dryfrac = VARIABLES.CANOPY.dryfrac;                                    %[] Canopy Dry fraction  
    wetfrac = VARIABLES.CANOPY.wetfrac;                                    %[] Canopy wet fraction  
    CAz = VARIABLES.CANOPY.CAz;                                            % [umol/mol] Atmospheric Concentration of CO2, for all layers 
    TAz = VARIABLES.CANOPY.TAz;                                            % [C] Air Temperature, for all layers 
    
    % VARIABLES    
    canstorheat = SWITCHES.canstorheat;                                    %[] Switch for considering canopy storage of heat [1] yes [0] no  
    
    % PARAMS
    ph_type = PARAMS.Photosyn.ph_type;                                     %[] Photosynthesis Type, 1 C3, 0 C4  
    
    % CONSTANTS
    Lv_kg = CONSTANTS.Lv_kg;                                               %[J / kg] Latent heat of vaporization                         

    % SWITCHES
    fsv_off = SWITCHES.fsv_off;
    
    % VERTSTRUC
    vinds_all = VERTSTRUC.vinds_all;                                       %[] Indicator of Layers where LAI is zero of lower for all species
    vinds = vinds_all{cntspecies};                                         %[] Indicator of Layers where LAI is zero of lower  
    
    
    % Initial Canopy States
    if (sunlit)  
%        Tl_prev =  VARIABLES.CANOPY.Tl_sun(:,cntspecies);
        Ci_prev =  VARIABLES.CANOPY.Ci_sun(:,cntspecies);                  % [umol/mol] Atmospheric Concentration of CO2 in the leaf in previous time step   
        gsv_prev =  VARIABLES.CANOPY.gsv_sun(:,cntspecies);                % [mol/m^2/s] Stomatal Conductance in previous time step    
    else
%        Tl_prev =  VARIABLES.CANOPY.Tl_shade(:,cntspecies);
        Ci_prev =  VARIABLES.CANOPY.Ci_shade(:,cntspecies);                % [umol/mol] Atmospheric Concentration of CO2 in the leaf in previous time step  
        gsv_prev =  VARIABLES.CANOPY.gsv_shade(:,cntspecies);              % [mol/m^2/s] Stomatal Conductance in previous time step 
    end
%*************************************************************************
%*************************************************************************
Tl = nan(size(dryfrac));
        

relax = 0;
relaxval = 0.25;
maxchange = 0.25;
maxiters = 10;
converged = 0;  cnt = 0;
while(~converged)
    
    % PHOTOSYNTHESIS
        if (ph_type(cntspecies)==1) 
            % C3
            [Ph, An, Ph_limit, Jc_C3, Jj_C3, Js_C3] = PHOTOSYNTHESIS_C3(VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, sunlit,cntspecies);
            
            % Set C4 Variables to NaN
            Jc_C4 = NaN(size(Ph)); Jj_C4 = NaN(size(Ph)); Js_C4 = NaN(size(Ph));
        else
            % C4
            [Ph, An ,Ph_limit, Jc_C4, Jj_C4, Js_C4] = PHOTOSYNTHESIS_C4(VARIABLES, PARAMS, VERTSTRUC, sunlit, cntspecies);
            
            % Set C3 Variables to NaN
            Jc_C3 = NaN(size(Ph)); Jj_C3 = NaN(size(Ph)); Js_C3 = NaN(size(Ph)); 
        end
        
        if (cnt>0)
            Ph = Ph - relax*(Ph-Ph_prev);
            An = An - relax*(An-An_prev);
        end
        
        if (sunlit)
            VARIABLES.CANOPY.An_sun(:,cntspecies) = An;
        else
            VARIABLES.CANOPY.An_shade(:,cntspecies) = An;
        end
    
    % LEAF BOUNDARY LAYER CONDUCTANCES
        [gbv, gbh] = BLC_Nikolov(VARIABLES, PARAMS, sunlit,cntspecies);
        if (sunlit)
            VARIABLES.CANOPY.gbv_sun(:,cntspecies) = gbv;
            VARIABLES.CANOPY.gbh_sun(:,cntspecies) = gbh;
        else
            VARIABLES.CANOPY.gbv_shade(:,cntspecies) = gbv;
            VARIABLES.CANOPY.gbh_shade(:,cntspecies) = gbh;
        end        
    
    % LEAF WATER POTENTIAL    
        [psil] = LEAF_WATER_POTENTIAL(VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, cntspecies);
        if (sunlit)
            VARIABLES.CANOPY.psil_sun(:,cntspecies) = psil;
        else
            VARIABLES.CANOPY.psil_shade(:,cntspecies) = psil;
        end
        [fsvg fsvm] = Tuzet_Function(VARIABLES, PARAMS, sunlit, cntspecies);
        if (fsv_off)
            if (sunlit)
                VARIABLES.CANOPY.fsvg_sun(:,cntspecies) = ones(size(fsvg));
                VARIABLES.CANOPY.fsvm_sun(:,cntspecies) = ones(size(fsvg));
            else
                VARIABLES.CANOPY.fsvg_shade(:,cntspecies) = ones(size(fsvm));
                VARIABLES.CANOPY.fsvm_shade(:,cntspecies) = ones(size(fsvm));                
            end
        else
            if (sunlit)
                VARIABLES.CANOPY.fsvg_sun(:,cntspecies) = fsvg;
                VARIABLES.CANOPY.fsvm_sun(:,cntspecies) = fsvm;
            else
                VARIABLES.CANOPY.fsvg_shade(:,cntspecies) = fsvg;
                VARIABLES.CANOPY.fsvm_shade(:,cntspecies) = fsvm;
            end
        end
     
    % STOMATAL CONDUCTANCE     
        [gsv, Ci, Hs, Cs] = BALL_BERRY(VARIABLES, PARAMS, VERTSTRUC, sunlit, cntspecies, SWITCHES, FORCING);
        if (cnt>0)
            gsv = gsv - relax*(gsv-gsv_prev);
            Ci = Ci - relax*(Ci-Ci_prev);
        end
        if (sunlit)
            VARIABLES.CANOPY.gsv_sun(:,cntspecies) = gsv;
            VARIABLES.CANOPY.Ci_sun(:,cntspecies) = Ci;
            VARIABLES.CANOPY.Cs_sun(:,cntspecies) = Cs;
            VARIABLES.CANOPY.Hs_sun(:,cntspecies) = Hs; 
        else
            VARIABLES.CANOPY.gsv_shade(:,cntspecies) = gsv;
            VARIABLES.CANOPY.Ci_shade(:,cntspecies) = Ci;
            VARIABLES.CANOPY.Cs_shade(:,cntspecies) = Cs;
            VARIABLES.CANOPY.Hs_shade(:,cntspecies) = Hs;
        end  

    % LEAF ENERGY BALANCE - DRY LEAF FRACTION
        if (canstorheat)
            [Tl_dry, H_dry, LE_dry, dHcan_dry, gv_dry] = LEB_QUARTIC_DRY_CHEAT (VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, sunlit,cntspecies);
        else
            [Tl_dry, H_dry, LE_dry, gv_dry] = LEB_QUARTIC_DRY (VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, sunlit,cntspecies);
            dHcan_dry = nan;            
        end
        
        if (sunlit)
                VARIABLES.CANOPY.TR_sun(:,cntspecies) = LE_dry/Lv_kg;    % [mm/s / LAI]  
            else
                VARIABLES.CANOPY.TR_shade(:,cntspecies) = LE_dry/Lv_kg;  % [mm/s / LAI]  
        end  
            

        TR = LE_dry / Lv_kg;    % [mm/s / LAI] 
                                                        
    
    % LEAF ENERGY BALANCE - WET LEAF FRACTION
    %   LE_wet = evaporation
    %   H_wet is alway zero
        if (canstorheat)    
            [Tl_wet, H_wet, LE_wet, dHcan_wet, gv_wet, vindswet, nvindswet] = LEB_QUARTIC_WET_CHEAT (VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, sunlit,cntspecies);
        else
            [Tl_wet, H_wet, LE_wet, gv_wet, vindswet, nvindswet] = LEB_QUARTIC_WET (VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, sunlit,cntspecies);
            dHcan_wet = nan;
        end            
            
    % Mean Leaf Temperature
        Tl(vindswet) = Tl_dry(vindswet).*dryfrac(vindswet) + Tl_wet(vindswet).*wetfrac(vindswet);
        Tl(nvindswet) = Tl_dry(nvindswet).*dryfrac(nvindswet);
    %    Tlabs = ((Tl_dry+273.15).^4.*dryfrac + (Tl_wet+273.15).^4.*wetfrac).^(1/4);
     %   Tl = Tlabs - 273.15;
        if (sunlit)
            VARIABLES.CANOPY.Tl_sun(:,cntspecies) = Tl;
        else
            VARIABLES.CANOPY.Tl_shade(:,cntspecies) = Tl;
        end  
        VARIABLES.CANOPY.Tl(:,cntspecies) = Tl;
    
        
    % TEST FOR SOLUTION DIVERGENCE
    if (cnt > 0)
        gsvdiff = gsv - gsv_prev;
        if isempty (vinds)
            break;
        else    
            % Check for solution divergence
            if ( max(abs((Ci(vinds)-Ci_prev(vinds))./Ci_prev(vinds))) > maxchange || ...
                 max(abs((gsv(vinds)-gsv_prev(vinds))./gsv_prev(vinds))) > maxchange || ...
                 max(abs((Tl(vinds)-Tl_prev(vinds))./Tl_prev(vinds))) > maxchange )

                % Rewind calculation and set relaxation on
                %gsv = gsv_prev;
                %Ci = Ci_prev;
                %Tl = Tl_prev;

                % Turn on relaxation
                relax = relaxval;

                gsv = gsv_prev - (1-relax)*(gsv_prev-gsv);
                Ci = Ci_prev - (1-relax)*(Ci_prev-Ci);
                Tl = Tl_prev - (1-relax)*(Tl_prev-Tl);



            elseif (relax > 0)
                relax = 0;
            end
        end        
        % Check for solution oscillation
        if (cnt>2)
           md1 = max(abs(gsvdiffprev));
           md2 = max(abs(gsvdiff));

           % Compare real numbers for "equality"
           if ( (md1 > (md2-0.01*md2)) && (md1 < (md2+0.01*md2)) )
               relax = relaxval;
               relax = 0;
           end
        end
        
        gsvdiffprev = gsvdiff;
    end
    
                                      
    % TEST CONVERGENCE
        if (cnt>0 & relax == 0)
            if ( (max(abs((gsv-gsv_prev)./gsv_prev)) < 0.01) && ...
                 (max(abs((Ci-Ci_prev)./Ci_prev)) < 0.01) && ...
                 (max(abs((Tl-Tl_prev)./Tl_prev)) < 0.01) )
                converged = 1;
            end
        end

        if (cnt>maxiters && converged==0)
            %disp(['*** TOO MANY INTERATIONS IN LEAF MODEL!!! --> Timestep:'])
            break;
        end

        % Update convergence check variables
        Ph_prev = Ph;
        An_prev = An;
        Ci_prev = Ci;
        gsv_prev = gsv;
        Tl_prev = Tl;
        
        cnt = cnt + 1;
end
% if cnt >= maxiters
%        Ph = real(Ph);
%        An = real(An);
%        Ci = real(Ci);
%        gsv = real(gsv);
%        Tl = real(Tl); 
%        H_dry = real(H_dry);
%        LE_dry = real(LE_dry);
%        LE_wet = real(LE_wet);
%        TR=real(TR);
% end       
% Compute energy fluxes, condensation and water storage for each layer
    H = H_dry.*dryfrac; % only dry fraction can produce sensible heat
    LE = LE_dry.*dryfrac + LE_wet.*wetfrac;
    
    dHcan = dHcan_dry.*dryfrac + dHcan_wet.*wetfrac;
    
% Compute Evaporation
    Evap_wm2 = LE_wet.*wetfrac;
    Evap_wm2(find(Evap_wm2<0)) = 0; % condensation instead of evaporation
    Evap_mm = Evap_wm2 / Lv_kg; % [mm/s / LAI]
  
% Compute Condensation
    Ch2o_mm_dry = LE_dry / Lv_kg; % [mm/s / LAI]
    Ch2o_mm_dry(find(Ch2o_mm_dry>0)) = 0;   % Transpiration instead of condensation
    
    Ch2o_mm_wet = LE_wet / Lv_kg; % [mm/s / LAI]
    Ch2o_mm_wet(find(Ch2o_mm_wet>0)) = 0;   % Evaporation instead of condensation
    
    Ch2o_mm = Ch2o_mm_dry + Ch2o_mm_wet; % [mm/s / LAI]

% Save number of iterations to reach convergence

VARIABLES.CANOPY.cntleafsolution(cntspecies) = cnt;
