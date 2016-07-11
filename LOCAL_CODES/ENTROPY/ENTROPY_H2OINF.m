function [SSresults] = ENTROPY_H2OINF (VARIABLES, VERTSTRUC, CONSTANTS, PARAMS, SSresults)

%=========================================================================
% This function computes the production of entropy that is associated  
% with the infiltration of water in the soil
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES       % VARIABLES structure
%       VERTSTRUC       % VARIABLES structure
%       CONSTANTS       % PARAMS structure
%       PARAMS          % CONSTANTS structure
%------------------------- Output Variables ------------------------------
%       SSresults       % VARIABLES structure
% 
%========================================================================


% Compute the entropy production due to wetting of water in the soil

% VARRIABLES 
smp = VARIABLES.SOIL.smp;                                                  %[mm] Soil Matric Potential
qlayer = VARIABLES.SOIL.qlayer;                                            %[mm/s] Fluxes of water between layers in the soil 
Ts = VARIABLES.SOIL.Ts;                                                    %[C] Temperature in the soil. For all layers  
hor = VARIABLES.SOIL.hor_drainage_lay;                                     %[mm/s] Horizontal Drainage 
qinfl = VARIABLES.SOIL.qinfl;                                              %[mm/s] Total Infiltration in the soil 
EVsoil = VARIABLES.SOIL.Esoil;                                              % Evaporation from soil [mm/s]

% VERTSTRUC
zmm=VERTSTRUC.znsmm;                                                       %[mm] Vector with depth at nodes of each layer

%CONSTANTS
mmH2OtoMPa = CONSTANTS.mmH2OtoMPa;                                         % Pressure conversion

%PARAMS
HC_liq = PARAMS.Soil.HC_liq;                                               % [J / kg / K] Specific Heat Capacity of liquid water 
% Fixing the constant soil layer problem
nl_soil=PARAMS.nl_soil;

% Change units
Hlayer = smp*mmH2OtoMPa - zmm*mmH2OtoMPa;                                  % Total Energy Head in [MPa] 
QQ = qlayer / 1000;                                                        % [m/s]  
INF =  (qinfl + EVsoil)/1000;                                              % [m/s] Infiltration into Soil      
HOR = hor/1000;                                                             % [m/s] Horizontal drainage  
BOT = QQ(end);                                                             % [m/s] Bottom Drainage 

TT = Ts + 273;                                                             % Soil Temperature [K] 

VHC_liq = HC_liq * 1000;                                                   % [J / m^3 / K] Volumetric heat capacity of liquid water 


% COMPUTE ENTROPY PRODUCTION WETTING 
% Fixing the constant soil layer problem
%for ii=2:1:12
for ii=2:1:nl_soil
    SH2O(ii-1) = ((Hlayer(ii-1) - Hlayer(ii))*QQ(ii))/((TT(ii-1)+TT(ii))/2)*10^6;                 % Entropy Production in [W/m^2/K]     
end

% COMPUTE ENERGY TOTAL THAT GOES INTO SOIL DUE TO WATER

% 1. Internal energy
%Eint = max(INF*VHC_liq*(TT(1)-273),0) - max(sum(HOR.*VHC_liq.*(TT-273)),0) ...
%    - max(BOT*VHC_liq*(TT(12)-273),0);                                     % Internal Energy Water [W/m2] 
% Fixing the constant soil layer problem
Eint = max(INF*VHC_liq*(TT(1)-273),0) - max(sum(HOR.*VHC_liq.*(TT-273)),0) ...
    - max(BOT*VHC_liq*(TT(nl_soil)-273),0);                                     % Internal Energy Water [W/m2] 


% Fixing the constant soil layer problem
%Ebond = (- BOT*Hlayer(12) - sum(HOR.*Hlayer))*10^6;                        % Bond Energy Water [W/m2]                        
Ebond = (- BOT*Hlayer(nl_soil) - sum(HOR.*Hlayer))*10^6;                        % Bond Energy Water [W/m2]                        

EH2O = Eint + Ebond;                                                       % Total Energy Budget [W/m2] 

SSresults.SH2O = SH2O;                                                     % [W/m2/K] Entropy due to water infiltration 
SSresults.EH2O = EH2O;                                                     % [W/m2] Energy associated with water budget in soil 
