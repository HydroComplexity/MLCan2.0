function [An_can, Ph_can, LE_can, H_can, dHcan, Rnrad_can, TR_can, ...
          Fc_soil, LE_soil, H_soil, Rnrad_soil, G, Tsurf, remainsoil, remaincan,remaineco, ...
          Rnrad_sun, Rnrad_shade, Rnrad_eco, ...
          An_sun, An_shade, LE_sun, LE_shade, H_sun, H_shade, TR_sun, TR_shade, ...
          Tl_sun, Tl_shade, psil_sun, psil_shade, gsv_sun, gsv_shade, fsvg_sun,fsvm_sun, ...
          fsvg_shade,fsvm_shade,Ci_sun, Ci_shade, CAz, TAz, EAz, Uz, gbv_sun, gbh_sun, gbv_shade, gbh_shade, ...
          LAIsun, LAIshade, fsun, fshade, ...
          Ph_limit_sun, Jc_C3_sun, Jj_C3_sun, Js_C3_sun, Jc_C4_sun, Jj_C4_sun, Js_C4_sun, ...
          Ph_limit_shade, Jc_C3_shade, Jj_C3_shade, Js_C3_shade, Jc_C4_shade, Jj_C4_shade, Js_C4_shade, ...
          PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, SWout, ...
          LWabs_can, LWemit_soil, LWemit_can, LWemit_sun, LWemit_shade, LWout, LWoutM, RH_soil, fdiff, ...
          Sh2o_prof, Sh2o_can, ppt_ground, Ch2o_prof, Ch2o_can, Evap_prof, Evap_can, ...
          dryfrac, wetfrac, Vz, VARIABLES, FORCING, ...
          SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out, SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
          LWabs_canM, LWabs_soilM, LSshaCON, LSsunCON] = ...    
        CANOPY_MODEL(SWITCHES, VERTSTRUC, FORCING, PARAMS, VARIABLES, CONSTANTS)
      
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % SWITCHES
    turb_on = SWITCHES.turb_on;                                            % Switch for simulation of turbulence in atmosphere [1] Yes [0] No     
    Sh2o_prof = VARIABLES.CANOPY.Sh2o_prof;                                % [mm] Intercepted water in canopy, for all layers     
    Tsoil = VARIABLES.SOIL.Ts(1);                                          % [C] Temperature in the Soil top layer 

    % VERTSTRUC
    vinds = VERTSTRUC.vinds;                                               % [] Indicator of Layers where LAI is zero of lower  
    nvinds = VERTSTRUC.nvinds;                                             % [] Indicator of Layers where LAI is zero of lower for all species         
    LAIz = VERTSTRUC.LAIz;                                                 % [m^2 leaf/ m^2 ground] LAI for all layers 
    LAIzall = VERTSTRUC.LAIzall;                                           % [m^2 leaf/ m^2 ground] LAI for all layers and species  
    fLAIz = VERTSTRUC.fLAIz;                                               % [] Fraction of total LAI for each species at each canopy layer 
    
    % PARAMS
    Ro = PARAMS.Resp.Ro;                                                   % [umol / m^2 / s] Ecosystem Respiration Parameter Q10 model
    Q10 = PARAMS.Resp.Q10;                                                 % [] Ecosystem respiration parameter Q10 model 
    nl_can = PARAMS.CanStruc.nl_can;                                       % [] Number of layers in canopy 
    nspecies = PARAMS.CanStruc.nspecies;                                   % [] Number of species 
    
    % CONSTANTS
    dtime = CONSTANTS.dtime;                                               % [s] Time step
    
    % Variables (Dongkook: For the output variables defined here)
    CAz = VARIABLES.CANOPY.CAz;                                            % [umol/mol] Atmospheric Concentration of CO2, for all layers
    TAz = VARIABLES.CANOPY.TAz;                                            % [C] Air Temperature, for all layers 
    EAz = VARIABLES.CANOPY.EAz;                                            % [kPa] Vapor Pressure Air, for all layers 
    
%*************************************************************************
%                          ALLOCATE MATRICES

%    ALLOCATE_STORAGE_CAN();
%percdiffprof = logical(zeros(nl_can,1));
%*************************************************************************
%*************************************************************************
             
% Vertical distribution of photosynthetic capacity
    [Vz] = PH_Dist(VERTSTRUC, PARAMS);
        VARIABLES.CANOPY.Vz = Vz;

% WIND PROFILE    
    [Uz, Km] = ORDER_1_CLOSURE_U(FORCING, VERTSTRUC, PARAMS);  
        VARIABLES.CANOPY.Uz = Uz;
        VARIABLES.CANOPY.Km = Km;
        
% INITIALIZE CANOPY STATES
    VARIABLES.CANOPY.gsv_sun = 0.01*ones(nl_can,nspecies);
    VARIABLES.CANOPY.gsv_shade = 0.01*ones(nl_can,nspecies);
        
    
% CANOPY PRECIPITATION INTERCEPTION
%   Smax [mm / LAI]
%   Sh2o_prof
%   ppt [mm]
    [Sh2o_prof, Smaxz, ppt_ground, wetfrac, dryfrac] = PRECIP_INTERCEPTION(FORCING, VARIABLES, VERTSTRUC, PARAMS);   
        VARIABLES.CANOPY.Sh2o_prof = Sh2o_prof;
        VARIABLES.CANOPY.Smaxz = Smaxz;
        VARIABLES.CANOPY.wetfrac = wetfrac;
        VARIABLES.CANOPY.dryfrac = dryfrac;
        VARIABLES.SOIL.ppt_ground = ppt_ground;

        % PREVIOUS
%TAz_prev = VARIABLES.CANOPY.TAz;
%VARIABLES_prev = VARIABLES;    
%diver_control = false;
    %======================================================================
    %                   SHORTWAVE RADIATION PROFILES
    %======================================================================
    [fsun, fshade, LAIsun, LAIshade, ...
     SWabs_sun, SWabs_shade, PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, ...
     PARabs_sun_lai, PARabs_shade_lai, ...
     SWabs_soil, PARabs_soil, NIRabs_soil, ...
     SWout, PARout, NIRout, PARtop, NIRtop, fdiff,...
     SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
     SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
     Kbm, taud, refl_soil] = ...
     SWRAD (FORCING, VERTSTRUC, PARAMS, CONSTANTS, VARIABLES); 
    
        % STORE VARIABLES
        VARIABLES.CANOPY.PARabs_sun = PARabs_sun_lai;
        VARIABLES.CANOPY.PARabs_shade = PARabs_shade_lai;
        VARIABLES.CANOPY.fsun = fsun;
        VARIABLES.CANOPY.fshade = fshade;
        VARIABLES.CANOPY.LAIsun = LAIsun;
        VARIABLES.CANOPY.LAIshade = LAIshade;
        VARIABLES.CANOPY.taud = taud;

    
    % LONGWAVE CONVERGENCE LOOP  
    converged_LW = 0; cnt_LW = 0; maxiters = 20; percdiff = 0.02;
    while (~converged_LW)     
        
        % LONGWAVE RADIATION ABSORPTION    
        [LWabs_can, LWabs_canM, LWabs_sun, LWabs_shade, LWabs_soil, LWabs_soilM, LWin, LWin2, LWout, LWoutM, ...
         LWemit_can, LWemit_sun, LWemit_shade, LWemit_soil] = ...
            LWRAD (FORCING, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);          

        % TOTAL ABSORBED RADIATION [W/m^2 ground]
        Totabs_sun = LWabs_sun + SWabs_sun;         
        Totabs_shade = LWabs_shade + SWabs_shade;   

        % TOTAL ABSORBED RADIATION PER UNIT LEAF AREA [W/m^2 leaf area]
        Rabs_sun_lai = Totabs_sun./LAIsun;
        Rabs_sun_lai(LAIsun==0)=0;
        Rabs_shade_lai = Totabs_shade./LAIshade;
        Rabs_shade_lai(LAIshade==0)=0;
        
        % SOIL ABSORBED ENERGY
        Totabs_soil = SWabs_soil + LWabs_soil;
        Rnrad_soil = Totabs_soil - LWemit_soil;

            % ASSIGN
            VARIABLES.CANOPY.Rnrad_soil = Rnrad_soil;
            VARIABLES.CANOPY.Rabs_sun = Rabs_sun_lai;
            VARIABLES.CANOPY.Rabs_shade = Rabs_shade_lai;
            VARIABLES.SOIL.Totabs_soil = Totabs_soil;
            FORCING.LWdn = LWin;
            FORCING.LWdn2 = LWin2;

        % ALLOCATE THE MEMORY FOR VARIABLES.    
        ALLOCATE_LEAF_VARS;    
        %==================================================================
        %                   SHADED CANOPY SOLUTION
        %   Calculations performed per [m^2 ground area], and canopy fluxes
        %   are calculated by integrating over the shaded leaf area
        %==================================================================     
        sunlit = 0;
        for ii=1:nspecies
         [Ph_shade(:,ii), An_shade(:,ii), Ci_shade(:,ii), gsv_shade(:,ii), Tl_shade(:,ii), LE_shade(:,ii),...
          dHcan_shade(:,ii), TR_shade(:,ii), Evap_shade(:,ii), H_shade(:,ii), ...
          psil_shade(:,ii), fsvg_shade(:,ii),fsvm_shade(:,ii), Ch2o_shade(:,ii), gbv_shade(:,ii), gbh_shade(:,ii), ...
          Ph_limit_shade(:,ii), Jc_C3_shade(:,ii), Jj_C3_shade(:,ii), Js_C3_shade(:,ii), Jc_C4_shade(:,ii),...
          Jj_C4_shade(:,ii), Js_C4_shade(:,ii), ...
          VARIABLES,LSshaCON(ii)] = ...
                     LEAF_SOLUTION(FORCING, VARIABLES, PARAMS, CONSTANTS, VERTSTRUC, SWITCHES, sunlit,ii);
        end 
        
       
                              
        %==================================================================
        %                       SUNLIT CANOPY SOLUTION
        %   Calculations performed per [m^2 leaf area], and canopy fluxes
        %       are calculated by integrating vertically over the sunlit 
        %       leaf area
        %==================================================================
        if (sum(fsun)==0)   % under nocturnal conditions all leaf area is 
                            % considered to be shaded
                            
            Ph_sun = zeros(nl_can,nspecies);
            An_sun = zeros(nl_can,nspecies);
            LE_sun = zeros(nl_can,nspecies);
            dHcan_sun = zeros(nl_can,nspecies);
            H_sun = zeros(nl_can,nspecies);
            Phtype_sun = NaN(nl_can,nspecies);
            gsv_sun = gsv_shade;
            gbv_sun = gbv_shade;
            gbh_sun = gbh_shade;
            Ci_sun = Ci_shade;
            VARIABLES.CANOPY.Cs_sun = zeros(nl_can,nspecies);
            Tl_sun = NaN(nl_can,nspecies);
            %Tl_sun = Tl_shade;   
            psil_sun = psil_shade;
            fsvg_sun = fsvg_shade;
            fsvm_sun = fsvm_shade;
            Evap_sun = zeros(nl_can,nspecies);
            TR_sun = zeros(nl_can,nspecies);
            Ch2o_sun = zeros(nl_can,nspecies);           
            
            Ph_limit_sun = zeros(nl_can,nspecies); 
            Jc_C3_sun = zeros(nl_can,nspecies); 
            Jj_C3_sun = zeros(nl_can,nspecies); 
            Js_C3_sun = zeros(nl_can,nspecies); 
            Jc_C4_sun = zeros(nl_can,nspecies); 
            Jj_C4_sun = zeros(nl_can,nspecies); 
            Js_C4_sun = zeros(nl_can,nspecies);
            LSsunCON = ones(1,nspecies);
            VARIABLES.CANOPY.Tl_sun = Tl_sun;
        else
            sunlit = 1;
            for ii=1:nspecies
            [Ph_sun(:,ii), An_sun(:,ii), Ci_sun(:,ii), gsv_sun(:,ii), Tl_sun(:,ii), LE_sun(:,ii),...
             dHcan_sun(:,ii), TR_sun(:,ii), Evap_sun(:,ii), H_sun(:,ii), ...
             psil_sun(:,ii), fsvg_sun(:,ii),fsvm_sun(:,ii), Ch2o_sun(:,ii), gbv_sun(:,ii), gbh_sun(:,ii), ...
             Ph_limit_sun(:,ii), Jc_C3_sun(:,ii), Jj_C3_sun(:,ii), Js_C3_sun(:,ii), Jc_C4_sun(:,ii),...
             Jj_C4_sun(:,ii), Js_C4_sun(:,ii),...
             VARIABLES,LSsunCON(ii)] = ...
                    LEAF_SOLUTION(FORCING, VARIABLES, PARAMS, CONSTANTS, VERTSTRUC, SWITCHES, sunlit,ii);                                
            end
        end
        % COMPUTE TEMPERATURE IN ALL CANOPY (AVERAGING OVER THE SPECIES)
            Tl_can_shade = nansum(((Tl_shade.*LAIzall)./repmat(LAIz,1,nspecies)),2);
            Tl_can_sun = nansum(((Tl_sun.*LAIzall)./repmat(LAIz,1,nspecies)),2);
        % CHECK THOSE LAYERS WHEN THERE IS NOT LAI AT ALL AND SET Tl_can = nan    
            Tl_can_shade(nvinds) = nan;
            Tl_can_sun(nvinds) = nan;
        % CHECK AGAIN THAT ALL THE SUNLIT FLUXES ARE nan IF fsun = 0
            inddark = fsun == 0;
            Tl_can_sun(inddark) = nan;
        % SAVE CANOPY TEMPERATURES
            VARIABLES.CANOPY.Tl_can_sun = Tl_can_sun;
            VARIABLES.CANOPY.Tl_can_shade = Tl_can_shade;
            
                            
        % ASSIGN
        
            VARIABLES.CANOPY.Ph_sun = Ph_sun;
            VARIABLES.CANOPY.Ph_shade = Ph_shade;        
            VARIABLES.CANOPY.An_sun = An_sun;
            VARIABLES.CANOPY.An_shade = An_shade;
            VARIABLES.CANOPY.LE_sun = LE_sun;
            VARIABLES.CANOPY.LE_shade = LE_shade;
            VARIABLES.CANOPY.H_sun = H_sun;
            VARIABLES.CANOPY.H_shade = H_shade;
        
            
        % SOIL RESPIRATION [umol CO2/ m^2 ground / s] 
            Fc_soil = Ro .* Q10.^((Tsoil - 10)/10); 
            
       % SOIL ENERGY FLUXES
       % determine the case
           zicesl = VARIABLES.SOIL.zicesl;         
         
        if SWITCHES.litter;         
           [H_soil, LE_soil, G, Gsl, RH_soil, Tsurf, remainsoil, dH, dS, VARIABLES] =...
           SOIL_SURFACE_FLUXES_LITTER(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES,FORCING);
        else        
           if zicesl > 0;                                  % Assumed ice to solve energy balance                               
                   [H_soil, LE_soil, G, Gsl, RH_soil, Tsurf, remainsoil, dH, dS, VARIABLES] =...
                   SOIL_SURFACE_FLUXES_SNOW(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES,FORCING);
           else                                            % Assumed no ice to solve energy balance
                  [H_soil, LE_soil, G, RH_soil, Tsurf, remainsoil, VARIABLES] = ... 
                  SOIL_SURFACE_FLUXES(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS,SWITCHES,FORCING);
                  dH = 0;
                  dS = 0;
           end
        end              
                % ASSIGN VARIABLES
                VARIABLES.SOIL.G = G;                
                VARIABLES.SOIL.Fc_soil = Fc_soil;
                VARIABLES.SOIL.LE_soil = LE_soil;
                VARIABLES.SOIL.H_soil = H_soil;
                VARIABLES.SOIL.Tsurf=Tsurf;
                VARIABLES.SOIL.remainsoil=remainsoil;
        
        if (turb_on)
            [CAz, EAz, TAz] = MICROENVIRONMENT(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);  
            % ASSIGN
                VARIABLES.CANOPY.CAz = CAz;
                VARIABLES.CANOPY.TAz = TAz;
                VARIABLES.CANOPY.EAz = EAz;
        end        
                       
                           
        % TEST LONGWAVE CONVERGENCE
        cnt_LW = cnt_LW + 1; 
        if (cnt_LW>1)
            diffprof = LWabs_can-LWabs_prev;
            percdiffprof(vinds) = diffprof(vinds)./LWabs_prev(vinds);
            
            Tl_new = [Tl_can_shade Tl_can_sun];
            difftl = max((max(abs(Tl_new-Tl_prev))));
            %difftl = 0;
            if  isempty(percdiffprof)
                converged_LW = 1;
                break
            end
            if (max(abs(percdiffprof(vinds))) < percdiff && difftl <= 0.5)
                converged_LW = 1;
            end
            if difftl > 10
                 break;
%                 VARIABLES = VARIABLES_prev;
%                 diver_control = true;
%                 cnt_LW = 0;
            end
%            if (diver_control && cnt_LW == 2)
%                break;
%            end
        end
        LWabs_prev = LWabs_can;
        Tl_prev = [Tl_can_shade Tl_can_sun];

        
        if (cnt_LW>maxiters && converged_LW==0)
            %disp(['*** TOO MANY ITERATIONS IN CANOPY MODEL!!! --> Timestep:', num2str(VARIABLES.niters_driver)]);
            break;
        end      
        
    end %LONGWAVE ITERATION
        
    
    % NET RADIATION
if (~isnan(FORCING.LWdn))
       Rnrad_eco = (FORCING.Rg - SWout) + (FORCING.LWdn - LWout);          %[W/m^2] Net Radiation in the Ecosystem 
else
       Rnrad_eco = (FORCING.Rg - SWout) + (LWin - LWout);                  %[W/m^2] Net Radiation in the Ecosystem                     
end
Rnrad_sun = SWabs_sun + LWabs_sun - LWemit_sun;                            %[W/m^2] Net Radiation Canopy Sunlit Fraction  
Rnrad_shade = SWabs_shade + LWabs_shade - LWemit_shade;                    %[W/m^2] Net Radiation Canopy Shade Fraction
        
    % H2O storage on foliage --> Precipitation and Condensation
Evap_prof = (sum(Evap_sun.*fLAIz,2).*LAIsun + sum(Evap_shade.*fLAIz,2).*LAIshade) * dtime;  % [mm] Evaporation at each layer
Evap_can = sum(Evap_prof);                                                 % [mm] Total canopy evaporation
        
Ch2o_prof = -(sum(Ch2o_sun.*fLAIz,2).*LAIsun + sum(Ch2o_shade.*fLAIz,2).*LAIshade) * dtime; % [mm] Condensation at each layer
Ch2o_can = sum(Ch2o_prof);                                                 % [mm] Total canopy condensation
   
   % ASSIGN
VARIABLES.CANOPY.Ch2o_prof = Ch2o_prof;                                    % [mm] Condensation Canopy at every Layer 
VARIABLES.CANOPY.Evap_prof = Evap_prof;                                    % [mm] Evaporation from Canopy at every Layer 
VARIABLES.CANOPY.Ch2o_can = Ch2o_can;                                      % [mm] Total Condensation Canopy 
VARIABLES.CANOPY.Evap_can = Evap_can;                                      % [mm] Total Evaporation from Canopy 
           
    % Adjust Canopy Water Storage
[Sh2o_prof, Sh2o_can, VARIABLES, Evap_can, Evap_prof, ppt_ground] = EVAP_CONDENSATION_ADJUST(VARIABLES, VERTSTRUC, PARAMS);
                        
    % COMPUTE CANOPY TOTAL FLUXES

Ph_can = sum(sum(Ph_sun.*fLAIz,2).*LAIsun) + sum(sum(Ph_shade.*fLAIz,2).*LAIshade);      %[umol/m^2/s] Total Photosynthesis canopy  
Ph_can_all = sum(Ph_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
            + sum(Ph_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[umol/m^2/s] Total Photosynthesis canopy for all  species
        
An_can = sum(sum(An_sun.*fLAIz,2).*LAIsun) + sum(sum(An_shade.*fLAIz,2).*LAIshade);      %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy 
An_can_all = sum(An_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
            + sum(An_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy for species
                
LE_can = sum(sum(LE_sun.*fLAIz,2).*LAIsun) + sum(sum(LE_shade.*fLAIz,2).*LAIshade);      %[W /m^2] Total Latent Heat canopy     
LE_can_all = sum(LE_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
            + sum(LE_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[W /m^2] Total Latent Heat canopy for all species  
        
LE_can_all_lay_shade = LE_shade.*fLAIz.*(repmat(LAIshade,1,nspecies));                   %[W /m^2] Total Latent Heat canopy for all layers and species, shade
LE_can_all_lay_sun = LE_sun.*fLAIz.*(repmat(LAIsun,1,nspecies));                         %[W /m^2] Total Latent Heat canopy for all layers and species, sunlit

H_can = sum(sum(H_sun.*fLAIz,2).*LAIsun) + sum(sum(H_shade.*fLAIz,2).*LAIshade);         %[W /m^2] Total Sensible Heat canopy   
H_can_all = sum(H_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
            + sum(H_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                        %[W /m^2] Total Sensible Heat canopy for all species     
        
H_can_all_lay_shade = H_shade.*fLAIz.*(repmat(LAIshade,1,nspecies));                     %[W /m^2] Total Sensible Heat canopy for all layers and species, shade
H_can_all_lay_sun = H_sun.*fLAIz.*(repmat(LAIsun,1,nspecies));                           %[W /m^2] Total Sensible Heat canopy for all layers and species, sunlit
        
        
dHcan = sum(sum(dHcan_sun.*fLAIz,2).*LAIsun) + sum(sum(dHcan_shade.*fLAIz,2).*LAIshade); %[W /m^2] Change in heat storage in canopy
dHcan_all = sum(dHcan_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...                       
            + sum(dHcan_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                    %[W /m^2] Change in heat storage in canopy for all species 
        
dHcan_all_lay_shade = dHcan_shade.*fLAIz.*(repmat(LAIshade,1,nspecies));                 %[W /m^2] Change in heat storage in canopy for all layers and species, shade fraction   
dHcan_all_lay_sun = dHcan_sun.*fLAIz.*(repmat(LAIsun,1,nspecies));                       %[W /m^2] Change in heat storage in canopy for all layers and species, sunlit fraction
                        
                               
TR_can = sum(sum(TR_sun.*fLAIz,2).*LAIsun) + sum(sum(TR_shade.*fLAIz,2).*LAIshade);      %[mm/s] = [g/m^2/s] Total Transpiration from Canopy
TR_can_all = sum(TR_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
            + sum(TR_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[mm/s] = [g/m^2/s] Total Transpiration from Canopy for all species
        
Rnrad_can = sum(Rnrad_sun(vinds)) + sum(Rnrad_shade(vinds));                             %[mm/s] Net Radiation Canopy 
        
     % COMPUTE REMAINDER
     
%        remaincan=Rnrad_eco - An_can*0.506 - H_can - LE_can - Totabs_soil;
remaincan = Rnrad_can - H_can - LE_can - An_can*0.506;                                   % [W / m^2] Energy Balance Error in Canopy                     
%        remaineco = Rnrad_eco - H_can - H_soil - LE_can - LE_soil - G - An_can*0.506 ;
remaineco = Rnrad_eco - (H_can + LE_can + An_can*0.506) ...
                        - (H_soil  + LE_soil + G + dH + dS);                             % [W / m^2] Energy Balance Error in all Ecosystem
%       remainsoil = Rnrad_soil - (H_soil + LE_soil + G + dH + dS);
        
if (Rnrad_eco-(Rnrad_can+Rnrad_soil)) > 3
     stop=34; 
end
        
%ASSIGN   
   
VARIABLES.CANOPY.Ph_can = Ph_can;                                          %[umol/m^2/s] Total Photosynthesis canopy
VARIABLES.CANOPY.Ph_can_all = Ph_can_all;                                  %[umol/m^2/s] Total Photosynthesis canopy for all  species   
        
VARIABLES.CANOPY.An_can = An_can;                                          %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy
VARIABLES.CANOPY.An_can_all = An_can_all;                                  %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy for species                          
        
VARIABLES.CANOPY.LE_can = LE_can;                                          %[W /m^2] Total Latent Heat canopy 
VARIABLES.CANOPY.LE_can_all = LE_can_all;                                  %[W /m^2] Total Latent Heat canopy for all species 
VARIABLES.CANOPY.LE_can_all_lay_shade = LE_can_all_lay_shade;              %[W /m^2] Total Latent Heat canopy for all layers and species, shade 
VARIABLES.CANOPY.LE_can_all_lay_sun = LE_can_all_lay_sun;                  %[W /m^2] Total Latent Heat canopy for all layers and species, sunlit 

        
VARIABLES.CANOPY.H_can = H_can;                                            %[W /m^2] Total Sensible Heat canopy 
VARIABLES.CANOPY.H_can_all = H_can_all;                                    %[W /m^2] Total Sensible Heat canopy for all species    
VARIABLES.CANOPY.H_can_all_lay_shade = H_can_all_lay_shade;                %[W /m^2] Total Sensible Heat canopy for all layers and species, shade 
VARIABLES.CANOPY.H_can_all_lay_sun = H_can_all_lay_sun;                    %[W /m^2] Total Sensible Heat canopy for all layers and species, sunlit 

VARIABLES.CANOPY.Tl_prev_dt = (nansum(Tl_sun.*fLAIz,2)).*fsun + ...
                                (nansum(Tl_shade.*fLAIz,2)).*fshade;       %[C] Temperature canopy for all layers  

VARIABLES.CANOPY.dHcan = dHcan;                                            %[W /m^2] Change in heat storage in canopy 
VARIABLES.CANOPY.dHcan_all = dHcan_all;                                    %[W /m^2] Change in heat storage in canopy for all species   
VARIABLES.CANOPY.dHcan_all_lay_shade = dHcan_all_lay_shade;                %[W /m^2] Total Sensible Heat canopy for all layers and species, shade 
VARIABLES.CANOPY.dHcan_all_lay_sun = dHcan_all_lay_sun;                    %[W /m^2] Total Sensible Heat canopy for all layers and species, sunlit 
                
VARIABLES.CANOPY.TR_can = TR_can;                                          %[mm/s] = [g/m^2/s] Total Transpiration from Canopy  
VARIABLES.CANOPY.TR_can_all = TR_can_all;                                  %[mm/s] = [g/m^2/s] Total Transpiration from Canopy for all species
VARIABLES.CANOPY.Evap_can = Evap_can;                                      %[mm] Total canopy evaporation
        
VARIABLES.CANOPY.Sh2o_prof = Sh2o_prof;                                    %[mm] Canopy moisture storage for all layers                                      
VARIABLES.CANOPY.Sh2o_can = Sh2o_can;                                      %[mm] Canopy canopy moisture storage  
        
VARIABLES.SOIL.Tlprev = VARIABLES.SOIL.Tli;                                % [C] Temperature Snow-Litter pack   
VARIABLES.refl_soil = refl_soil;                                           % [] SW Reflection from Soil
        
                
                 