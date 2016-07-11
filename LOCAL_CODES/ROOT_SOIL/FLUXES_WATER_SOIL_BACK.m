function [VARIABLES] = FLUXES_WATER_SOIL_BACK (VARIABLES, VERTSTRUC, PARAMS, CONSTANTS)

%=========================================================================
% This functions corrects the fluxes of water in the surface. 
% The implicit solution may fail when the infiltration rate is too high
% This function corrects for that and allocates the fluxes into the litter
% layer or snow pack or into runoff
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES       % VARIABLES structure
%       VERTSTRUC       % VARIABLES structure
%       PARAMS          % PARAMS structure
%       CONSTANTS       % CONSTANTS structure
%------------------------- Output Variables ------------------------------
%       VARIABLES       % VARIABLES structure
% 
%========================================================================


% Dereference blocks

% CONSTANTS
dtime = CONSTANTS.dtime;                                                    % [s] Time Step

% VERTSTRUC
porsl = VERTSTRUC.porsl;                                                    % [] Porosity
dzsmm = VERTSTRUC.dzsmm;                                                    % [mm] Soil Layer Tickness  

% PARAMS
rho_liq = PARAMS.Soil.rho_liq;                                              % [kg / m3] Water Density
HC_snow = PARAMS.Soil.HC_snow;                                              % [] Maximum fraction water capacity that snow can hold 

% VARIABLES
dzlit_mm = VARIABLES.SOIL.dzlit_m*1000;                                     % [mm] Thickness of litter 
qinflL = VARIABLES.SOIL.qinflL;                                             % [mm/s] Water that infiltrates into snow
net_qinflL = VARIABLES.SOIL.net_qinflL;                                     % [mm/s] Net water that infiltrates into snow excluding drainage in soil
volliq = VARIABLES.SOIL.volliq;                                             % [] Volumetric water content
qinfl = VARIABLES.SOIL.qinfl;                                               % [mm/s] Net Infiltration into the soil
qlayer = VARIABLES.SOIL.qlayer;                                             % [mm/s] Fluxes of water between layers in the soil
wliqsl = VARIABLES.SOIL.wliqsl;                                             % [kg / m^2] Liquid Water Density per unit area
wicesl = VARIABLES.SOIL.wicesl;                                             % [kg / m^2] Ice Water Density per unit area
wsn = VARIABLES.SOIL.wsn;                                                   % [kg / m^2] Snow Water Density per unit area
zliqsl = VARIABLES.SOIL.zliqsl;                                             % [mm] Liquid water depth
zicesl = VARIABLES.SOIL.zicesl;                                             % [mm] Ice water depth
zsn = VARIABLES.SOIL.zsn;                                                   % [mm] Snow depth
dwat = VARIABLES.SOIL.dwat;                                                % [] Change in Soil Moisture
nstype = VARIABLES.SOIL.type;                                              % Type os solution implemented in the soil


% CCOMPUTE QBACK 
% This will occur if  
%1)Richards equation can not be solved due to a high flux top BC
%2)If flux top BC is higher than Ksa
qback = max(qinfl - qlayer(1),0);                                                  %[mm/s] Flux of water returned to litter or surface by sypersaturation of first layer
%However, the return to soil to saturate top layer occurs only in condition 1) 
if nstype(2) == 1
    qadflux = min((porsl(1)-volliq(1))*dzsmm(1)/dtime,qback);
else
    qadflux = 0;
end
    
volliq(1) = min(volliq(1) + qadflux*dtime/dzsmm(1) , porsl(1)); 
dwat(1) = min(dwat(1) + qadflux*dtime/dzsmm(1),porsl(1)); 

qinfl = qlayer(1) + qadflux;
zback_res = (qback - qadflux)*dtime;
       
% Compute how much of the flux back could be store in the litter or
% snowpack
if zsn > 0                                                                 % [mm]
    dz = min(max((HC_snow-wliqsl/wsn)*zsn,0),zback_res);
    zliqsl_new = dz + zliqsl;
    volliqli = 1;                                                          % [-]                    
    % check (qback = qadflux + qback_lit + runoff)
    qback_lit = dz/dtime;                                                  % Compute the flux that goes back to litter 
    runoff = qback - qadflux - qback_lit;                                  % Compute Runoff Flux [mm/s]         

else 
    zliqsl_new = 0;
    qback_lit = 0;
    runoff = qback - qadflux + (zliqsl)/dtime;
    volliqli = 0;
end
        
qinflL = qinflL - runoff;
net_qinflL = net_qinflL + qback_lit;
        
% Update states
wliqsl_new = (zliqsl_new/1000)*rho_liq;                                    % [kg/m2]
wsn = wicesl + wliqsl_new;                                                 % [kg/m2] 
        
% ASSIGN VALUES FOR LITTER
VARIABLES.SOIL.volliq = volliq;                                            % [] Soil Volumetric water content 
VARIABLES.SOIL.volliqli = volliqli;                                        % [] Litter Volumetric water content 
VARIABLES.SOIL.qinfl = qinfl;                                              % [mm/s] Net Infiltration into the soil 
VARIABLES.SOIL.dwat = dwat;                                                % [] Change in volumetric water content 

VARIABLES.SOIL.qback = qback;                                              % [mm/s] Flux returning to snow-litter pack 
VARIABLES.SOIL.qadflux = qadflux;                                          % [mm/s] Final Flux that Infiltrates in the soil 
VARIABLES.SOIL.qback_lit = qback_lit;                                      % [mm/s] Real flux that is going back to snow 
VARIABLES.SOIL.runoff = runoff;                                            % [mm/s] Runnof Flux 

VARIABLES.SOIL.qinflL = qinflL;                                            % [mm/s] Infiltration Flux Litter 
VARIABLES.SOIL.net_qinflL = net_qinflL;                                    % [mm/s] Net Infiltration Flux  
        
VARIABLES.SOIL.zliqsl = zliqsl_new;                                        % [mm] Liquid water depth  
VARIABLES.SOIL.wliqsl = wliqsl_new;                                        % [kg / m^2] Liquid Water Density per unit area 
VARIABLES.SOIL.wsn = wsn;                                                  % [kg / m^2] Snow Water Density per unit area        
