function [VARIABLES] = MASS_BALANCE (VARIABLES, CONSTANTS, PARAMS, FORCING, SWITCHES, VERTSTRUC, tt) 


%=========================================================================
% THIS FUNCTION HELPS TO CHECK THE WATER BALANCE IN THE CODE.
% IT COMPUTES THE WATER MASS BALANCE IN: 
% i.   Canopy
% ii.  The soil
% iii. The litter-soil 
% iv.  The litter-soil-canopy
% THE RESULTS ARE SAVE IN THE STRUCTURE VARIABLES.MB
% Dereference bloxk
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES       % SWITCHES structure
%       CONSTANTS       % CONSTANTS structure
%       PARAMS          % PARAMS structure
%       FORCING         % FORCING structure
%       SWITCHES        % SWITCHES structure
%       VERTSTRUC       % VERTSTRUC structure
%       tt              % Time Step
%------------------------- Output Variables ------------------------------
%       VARIABLES       % SWITCHES structure
% 
%========================================================================
% Dereference blocks

% CONSTANTS
dtime = CONSTANTS.dtime;                                                    % [s] Delta Time  

% FORCING
PPTf = FORCING.ppt/dtime;                                                   % [mm/s] Total Rainfall  

% VARIABLES
dwat = VARIABLES.SOIL.dwat;                                                 % [] Change in Soil Moisture
qinfl = VARIABLES.SOIL.qinfl;                                               % [mm/s] Total Infiltration in the soil 
PPTgr = VARIABLES.SOIL.pptrate_ground;                                      % [mm/s] Total Rainfall Reaching Ground  
EVsoil = VARIABLES.SOIL.Esoil;                                              % [mm/s] Evaporation from soil 
Esl = VARIABLES.SOIL.Esl;                                                   % [mm/s] Evaporation from soil 
EVcanf = VARIABLES.CANOPY.Evap_can/dtime;                                   % [mm/s] Evaporation from canopy 
qlayer = VARIABLES.SOIL.qlayer;                                             % [mm/s] Fluxes in the soil 
Bdrainage = qlayer(end);                                                    % [mm/s] Bottom Drainage from soil 
Hdrainage = VARIABLES.SOIL.hor_drainage;                                    % [mm/s] Horizontal Drainage from soil      

% Dongkook
flux_Ss = VARIABLES.SOIL.flux_Ss;                                           % [1/s] Flux Moisture due to Storativity
% Dongkook - End

TR = VARIABLES.CANOPY.TR_can_all;                                           % [mm/s] Transpiration all species 
runoff = VARIABLES.SOIL.runoff;                                             % [mm/s] Runoff species 
Sh2o_can = VARIABLES.CANOPY.Sh2o_can;                                       % [mm] Total water content in the canopy 
Sh2o_can_prev = VARIABLES.CANOPY.Sh2o_can_prev;                             % [mm] Total water content in the canopy in previous time spte 
COcanf = VARIABLES.CANOPY.Ch2o_can/dtime;                                   % [mm/s] Total canopy condensation 
zliqsl = VARIABLES.SOIL.zliqsl;                                             % [mm] Liquid water depth  in Snow-litter pack     
zicesl = VARIABLES.SOIL.zicesl;                                             % [mm] Ice depth  in Snow-litter pack      
zliqsl_prev = VARIABLES.SOIL.zliqsl_prev;                                   % [mm] Liquid water depth  in Snow-litter pack     
zicesl_prev = VARIABLES.SOIL.zicesl_prev;                                   % [mm] Ice depth  in Snow-litter pack     
zliq_gro = VARIABLES.SOIL.zliq_gro;                                         % [mm] Liquid water reaching snow-litter pack
zice_gro = VARIABLES.SOIL.zice_gro;                                         % [mm] Ice water reaching snow-litter pack
Esl_ice = VARIABLES.SOIL.Esl_ice;                                           % [mm/s] Sublimation from snow-litter pack
Esl_liq = VARIABLES.SOIL.Esl_liq;                                           % [mm/s] Evaporation from snow-litter pack

% PARAMS
rho_liq = PARAMS.Soil.rho_liq;                                              % [kg/m3] Density water
rho_ice = PARAMS.Soil.rho_ice;                                              % [kg/m3] Density Ice

% VERTSTRUC
dzsmm = VERTSTRUC.dzsmm;                                                    % [mm] Vector with Layer thicknesses


%INFsoil =  drainlitter/dtime;                                              % [mm/s] Infiltration into Soil, Drainage From Litter 
INFsoil =  qinfl + EVsoil;                                                  % [mm/s] Infiltration into Soil       
qadflux = VARIABLES.SOIL.qadflux;                                           % Additional flow that is going to soil by drainage
                                                                            % in those scenarions when the implicit solution can not
                                                                            % compute the infiltration correctly 

% change in states
if tt>1
    dws_f = sum(dwat.*dzsmm)/dtime;                                         % [mm/s] Rate of change in soil water content 
    % Dongkook
    dwSs_f = sum(flux_Ss.*dzsmm);                                           % [mm/s] Rate of change in soil water content
    % Dongkook - End    

    dwc_f = (Sh2o_can - Sh2o_can_prev)/dtime;                               % [mm/s]Rate of change in canopy water content 
    

    % COMPUTE CANOPY MASS BALANCE
    inputmb = PPTf + COcanf;
    outputmb = EVcanf + PPTgr;
    dstate = dwc_f;

    mbcan = (inputmb - outputmb - dstate)*86400; %[mm/day]    
        
    if SWITCHES.litter
        % COMPUTE SOIL MASS BALANCE
        inputmb = INFsoil;
        outputmb = sum(TR) + EVsoil + Bdrainage + Hdrainage ;% INFqback;
        % Dongkook
        %dstate = dws_f;  
        dstate = dws_f + dwSs_f;   
        % Dongkook - End
        mbsoil = (inputmb - outputmb - dstate)*86400; %[mm/day]   
                        
        % COMPUTE SOIL - LITTER SNOW PACK MASS BALANCE
        dwl_f = (zliqsl - zliqsl_prev + (zicesl - zicesl_prev)*(rho_ice/rho_liq))/dtime;
        inputmb = (zliq_gro + zice_gro*(rho_ice/rho_liq))/dtime;
        outputmb = sum(TR) + EVsoil + Esl_ice*(rho_ice/rho_liq) + Esl_liq + Bdrainage + Hdrainage + runoff;
        % Dongkook
        %dstate = dws_f + dwl_f;       
        dstate = dws_f + dwl_f + dwSs_f; 
        % Dongkook - End
        mblittersoil = (inputmb - outputmb - dstate)*86400; %[mm /day]

        % COMPUTE CANOPY-SOIL-LITTER MASS BALANCE
        inputmb = PPTf + COcanf;
        outputmb = EVcanf + sum(TR) + EVsoil + Esl_ice*(rho_ice/rho_liq) + Esl_liq + Bdrainage + Hdrainage + runoff;
        % Dongkook
        %dstate = dws_f + dwl_f + dwc_f;   
        dstate = dws_f + dwl_f + dwc_f + dwSs_f; 
        % Dongkook - End
        mbcanlittersoil = (inputmb - outputmb - dstate)*86400; %[mm/day] 
    else
        % COMPUTE SOIL MASS BALANCE
        inputmb = INFsoil;
        outputmb = sum(TR) + EVsoil + Bdrainage + Hdrainage ;% INFqback;
        % Dongkook
        %dstate = dws_f;  
        dstate = dws_f + dwSs_f;   
        % Dongkook - End
        mbsoil = (inputmb - outputmb - dstate)*86400; %[mm/day]   
                        
        % COMPUTE SOIL - LITTER SNOW PACK MASS BALANCE
        dwl_f = (zliqsl - zliqsl_prev + (zicesl - zicesl_prev)*(rho_ice/rho_liq))/dtime;
        inputmb = (zliq_gro + zice_gro*rho_ice/rho_liq)/dtime;
        outputmb = sum(TR) + EVsoil + Esl_ice*(rho_ice/rho_liq) + Bdrainage + Hdrainage + runoff;
        % Dongkook
        %dstate = dws_f + dwl_f;  
        dstate = dws_f + dwl_f + dwSs_f;  
        % Dongkook - End
        mblittersoil = (inputmb - outputmb - dstate)*86400; %[mm/day]

        % COMPUTE CANOPY-SOIL-LITTER MASS BALANCE
        inputmb = PPTf + COcanf;
        outputmb = EVcanf + sum(TR) + EVsoil + Esl_ice*(rho_ice/rho_liq) + Bdrainage + Hdrainage + runoff;
        % Dongkook
        %dstate = dws_f + dwl_f + dwc_f; 
        dstate = dws_f + dwl_f + dwc_f + dwSs_f;  
        % Dongkook - End
        mbcanlittersoil = (inputmb - outputmb - dstate)*86400; %[mm/day] 
    end

% SAVE VARIABLES IN STRUCTURE
% VARIABLES.MB.PPTf = PPTf;                                                  %[mm/s] Total Rainfall  
% VARIABLES.MB.PPTgr = PPTgr;                                                %[mm/s] Total Rainfall Reaching Ground   
% VARIABLES.MB.INFsoil = INFsoil;                                            %[mm/s] Infiltration into Soil
% VARIABLES.MB.EVsoil = EVsoil;                                              %[mm/s] Evaporation from soil 
% VARIABLES.MB.EVcanf = EVcanf;                                              %[mm/s] Evaporation from canopy 
% VARIABLES.MB.Bdrainage = Bdrainage;                                        %[mm/s] Bottom Drainage from soil
% VARIABLES.MB.Hdrainage = Hdrainage;                                        %[mm/s] Horizontal Drainage from soil
% VARIABLES.MB.TR = TR;                                                      %[mm/s] Transpiration all species 
% VARIABLES.MB.runoff = runoff;                                              %[mm/s] Runoff 
% VARIABLES.MB.COcanf = COcanf;                                              %[mm/s] Total canopy condensation 
% VARIABLES.MB.dws_f = dws_f;                                                %[mm/s] Rate of change in soil water content  
% VARIABLES.MB.dwc_f = dwc_f;                                                %[mm/s] Rate of change in canopy water content 
% VARIABLES.MB.dwl_f = dwl_f;                                                %[mm/s] Rate of change in snow-litter water content 
VARIABLES.MB.mbsoil = mbsoil;                                              %[mm/day] Mass Balance Error Soil                                
VARIABLES.MB.mbcan = mbcan;                                                %[mm/day] Mass Balance Error Canopy
VARIABLES.MB.mblittersoil = mblittersoil;                                  %[mm/day] Mass Balance Error Litter and Soil 
VARIABLES.MB.mbcanlittersoil = mbcanlittersoil;                            %[mm/day] Mass Balance Error Canopy-Litter-Soil 

end



