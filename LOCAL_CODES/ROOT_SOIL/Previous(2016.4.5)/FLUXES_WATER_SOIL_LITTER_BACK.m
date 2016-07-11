function [VARIABLES] = FLUXES_WATER_SOIL_LITTER_BACK (VARIABLES, VERTSTRUC, PARAMS, CONSTANTS)

% This functions corrects the fluxes of water in the surface. 
% The implicit solution may fail when the infiltration rate is too high
% This function corrects for that and allocates the fluxes into the litter
% layer or snow pack or into runoff


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
% CONSTANS
dtime = CONSTANTS.dtime;                                                    % [s] Time Step

% VERTSTRUC
porsl = VERTSTRUC.porsl;                                                    % Soil Porosity
dzsmm = VERTSTRUC.dzsmm;                                                    % [mm] Soil Layer Tickness 

% PARAMS
thetals = PARAMS.Soil.thetals;                                             % [] Volumetric Water Content in Litter at Saturation 
rho_liq = PARAMS.Soil.rho_liq;                                             % [kg / m^3] Density Liquid Water 
rho_ice = PARAMS.Soil.rho_ice;                                             % [kg / m^3] Density Ice Water 
HC_snow = PARAMS.Soil.HC_snow;                                             % [] Maximum fraction water capacity that snow can hold 

% VARIABLES
case1sl = VARIABLES.SOIL.case1sl;                                          % Index that defines whether is treated as snow or as a litter
dzlit_mm = VARIABLES.SOIL.dzlit_m*1000;                                    % Thickness of litter [mm]
qinflL = VARIABLES.SOIL.qinflL;                                            % [mm/s] Water that infiltrates into snow 
net_qinflL = VARIABLES.SOIL.net_qinflL;                                    % [mm/s] 
volliq = VARIABLES.SOIL.volliq;                                            % Volumetric Water Content  
qinfl = VARIABLES.SOIL.qinfl;                                              % [mm/s] Net Infiltration into the soil 
qlayer = VARIABLES.SOIL.qlayer;                                            % [mm/s] Fluxes of water between layers in the soil 
wliqsl = VARIABLES.SOIL.wliqsl;                                            % [kg / m^2] Liquid Water Density per unit area 
wicesl = VARIABLES.SOIL.wicesl;                                            % [kg / m^2] Ice Water Density per unit area 
wsn = VARIABLES.SOIL.wsn;                                                  % [kg / m^2] Snow Water Density per unit area 
zliqsl = VARIABLES.SOIL.zliqsl;                                            % [mm] Liquid water depth 
zicesl = VARIABLES.SOIL.zicesl;                                            % [mm] Ice water depth 
zsn = VARIABLES.SOIL.zsn;                                                  % [mm] Snow depth 


% CCOMPUTE QBACK 
% This can IF 
%1)Richards equation can not be solved due to a high flux top BC
%2)If flux top BC is higher than Ksa

%qback = qinfl - qlayer(1);
%qadflux = 0;
%qadflux = min((porsl(1)-volliq(1))*dzsmm(1)/dtime,qback);
%volliq(1) = porsl(1);
%qinfl = qlayer(1) + qadflux;
%zback_res = (qback - qadflux)*dtime;
       


qback = qinfl - qlayer(1);                                                  %[mm/s] Flux of water returned to litter or surface by sypersaturation of first layer
    
qadflux = min((porsl(1)-volliq(1))*dzsmm(1)/dtime,qback);
volliq(1) = min(volliq(1) + qadflux*dtime/dzsmm(1) , porsl(1)); 
qinfl = qlayer(1) + qadflux;
zback_res = (qback - qadflux)*dtime;

% Compute how much of the flux back could be store in the litter or
% snowpack
if (zliqsl < dzlit_mm*thetals && case1sl==0 )                              % No Ice at all  
    dz = max(min(zback_res, dzlit_mm*thetals  - zliqsl - zicesl),0);
    zliqsl_new = dz + zliqsl;
elseif (zliqsl > dzlit_mm*thetals && case1sl==0)                           % No ice oversaturated 
    dz = 0;
    %zliqsl_new = dz + zliqsl;
    zliqsl_new = dzlit_mm*thetals;    
    qback = qback + (zliqsl - dzlit_mm*thetals)/dtime;                     % all the excess over the litter goes to runoff
elseif (zsn < dzlit_mm*thetals && case1sl==1)                              % Ice but snowpack is lower than litter 
    dz = min(max(dzlit_mm*thetals - zliqsl - zicesl,0)+max((HC_snow-wliqsl/wsn)*zsn,0),zback_res);
    zliqsl_new = dz + zliqsl;
elseif (zsn > dzlit_mm*thetals && case1sl==1)                              % Ice and snowpack is higher than litter 
    dz = min(max((HC_snow-wliqsl/wsn)*zsn,0),zback_res);
    zliqsl_new = dz + zliqsl;
end
volliqli = zliqsl_new/dzlit_mm;                                            % [-]                    
if volliqli < 0
    volliqli = PARAMS.Soil.thetamin;
end
% check (qback = qadflux + qback_lit + runoff)
qback_lit = dz/dtime;                                                      % Compute the flux that goes back to litter 
runoff = qback - qadflux - qback_lit;                                      % Compute Runoff Flux [mm/s]         
        
qinflL = qinflL - runoff;
net_qinflL = net_qinflL + qback_lit;
        
% Update states
wliqsl = (zliqsl_new/1000)*rho_liq;                                        %[kg/m2]
wsn = wicesl + wliqsl;                                                     % [kg/m2] 
        
% ASSIGN VALUES FOR LITTER
VARIABLES.SOIL.volliq = volliq;                                            % [] Soil Volumetric water content 
VARIABLES.SOIL.volliqli = volliqli;                                        % [] Litter Volumetric water content 
VARIABLES.SOIL.qinfl = qinfl;                                              % [mm/s] Net Infiltration into the soil 

VARIABLES.SOIL.qback = qback;                                              % [mm/s] Flux returning to snow-litter pack 
VARIABLES.SOIL.qadflux = qadflux;                                          % [mm/s] Final Flux that Infiltrates in the soil 
VARIABLES.SOIL.qback_lit = qback_lit ;                                     % [mm/s] Real flux that is going back to snow-liter pack. 
VARIABLES.SOIL.runoff = runoff;                                            % [mm/s] Runnof Flux                    

VARIABLES.SOIL.qinflL = qinflL;                                            % [mm/s] Infiltration Flux Litter 
VARIABLES.SOIL.net_qinflL = net_qinflL;                                    % [mm/s] Net Infiltration Flux     
        
VARIABLES.SOIL.zliqsl = zliqsl_new;                                        % [mm] Liquid water depth 
VARIABLES.SOIL.wliqsl = wliqsl;                                            % [kg / m^2] Liquid Water Density per unit area 
VARIABLES.SOIL.wsn = wsn;                                                  % [kg / m^2] Snow Water Density per unit area 
