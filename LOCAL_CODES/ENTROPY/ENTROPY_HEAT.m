function [SSresults] = ENTROPY_HEAT (VARIABLES, SWITCHES, SSresults)

%=========================================================================
% This code computes the fluxes of entropy due to heat. 
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES                      % VARIABLES structure
%       SWITCHES                       % SWITCHES structure
%       SSresults                      % SSresults structure 
%------------------------- Output Variables ------------------------------
%       SSresults                      % SSresults structure 
% 
%========================================================================
% De reference block
%*********************************************************************

LAIsun = VARIABLES.CANOPY.LAIsun;                                          %[m^2 leaf / m^2 ground] LAI in sunlit   
LAIshade = VARIABLES.CANOPY.LAIshade;                                      %[m^2 leaf / m^2 ground] LAI in shade 

LE_soil = VARIABLES.SOIL.LE_soil;                                          %[W/m^2] Latent Heat from soil 
H_soil = VARIABLES.SOIL.H_soil;                                            %[W/m^2] Sensible Heat from soil 
dH = VARIABLES.SOIL.dH;                                                    %[W/m2] Delta of Energy Melting/Fusion 
dS = VARIABLES.SOIL.dS;                                                    %[W/m2] Change of Internal Energy in Snow-Litter Pack in Time Step 

if SWITCHES.litter;
    G = VARIABLES.SOIL.G;                                                  %[W(m^2] Ground Heat Flux  
else
    G = VARIABLES.SOIL.G;                                                  %[W(m^2] Ground Heat Flux 
end
plants = SWITCHES.plants;

LE_can_all_lay_shade = VARIABLES.CANOPY.LE_can_all_lay_shade;              %[W/m^2] LE from shade canopy for all layers and for all species
LE_can_all_lay_sun = VARIABLES.CANOPY.LE_can_all_lay_sun;                  %[W/m^2] LE from sunlit canopy for all layers and for all species 
LE_can_all_lay = LE_can_all_lay_shade + LE_can_all_lay_sun;                %[W/m^2] LE from canopy for all layers and for all species 

H_can_all_lay_shade = VARIABLES.CANOPY.H_can_all_lay_shade;                %[W/m^2] H from shade canopy for all layers and for all species 
H_can_all_lay_sun = VARIABLES.CANOPY.H_can_all_lay_sun;                    %[W/m^2] H from sunlit canopy for all layers and for all species 
H_can_all_lay = H_can_all_lay_shade + H_can_all_lay_sun;                   %[W/m^2] H from canopy for all layers and for all species 

dHcan_all_lay_shade = VARIABLES.CANOPY.dHcan_all_lay_shade;                %[W/m^2] Change in internal heat content in shade leaves, for all layers and species  
dHcan_all_lay_sun = VARIABLES.CANOPY.dHcan_all_lay_sun;                    %[W/m^2] Change in internal heat content in sunlit leaves, for all layers and species 
dHcan_all_lay = dHcan_all_lay_shade + dHcan_all_lay_sun;                   %[W/m^2] Change in internal heat content in leaves, for all layers and specie 


% change temperature units from C to K
Tsurf = VARIABLES.SOIL.Tsurf + 273.15;                                     %[K] Soil Surface Temperature 
Tl_sun = VARIABLES.CANOPY.Tl_sun + 273.15;                                 %[K] Sunlit Leaf Temperature for all layers and species
Tl_shade = VARIABLES.CANOPY.Tl_shade + 273.15;                             %[K] Shade Leaf Temperature for all layers and species
Tl_can_sun = VARIABLES.CANOPY.Tl_can_sun + 273.15;                         %[K] Sunlit Leaf Temperature for all layers, averaged through species                          
Tl_can_shade = VARIABLES.CANOPY.Tl_can_shade + 273.15;                     %[K] Shade Leaf Temperature for all layers, averaged through species

Tl_net = SSresults.Tl_net;                                                 %[K] Net (All Ecosystem) Temperature 1 Based on LW 


%**************************************************************************
%**************************************************************************
%                         I.  FOR THE CANOPY
%**************************************************************************
%**************************************************************************

% Latent Heat
SSLE_can_all_lay_shade = LE_can_all_lay_shade./Tl_shade;
SSLE_can_all_lay_shade(isnan(SSLE_can_all_lay_shade)) = 0;
SSLE_can_all_lay_sun = LE_can_all_lay_sun./Tl_sun;
SSLE_can_all_lay_sun(isnan(SSLE_can_all_lay_sun)) = 0;


SSLE_can_all_lay = SSLE_can_all_lay_shade + SSLE_can_all_lay_sun;
SSLE_can_all = sum(SSLE_can_all_lay);        
SSLE_can = sum(SSLE_can_all);        

% Sensible Heat
SSH_can_all_lay_shade = H_can_all_lay_shade./Tl_shade;
SSH_can_all_lay_shade(isnan(SSH_can_all_lay_shade)) = 0;
SSH_can_all_lay_sun = H_can_all_lay_sun./Tl_sun;
SSH_can_all_lay_sun(isnan(SSH_can_all_lay_sun)) = 0;

SSH_can_all_lay = SSH_can_all_lay_shade + SSH_can_all_lay_sun;
SSH_can_all = sum(SSH_can_all_lay);        
SSH_can = sum(SSH_can_all);        

% Leaf Heat Storage
SSdHcan_all_lay_shade = dHcan_all_lay_shade./Tl_shade;
SSdHcan_all_lay_shade(isnan(SSdHcan_all_lay_shade)) = 0;
SSdHcan_all_lay_sun = dHcan_all_lay_sun./Tl_sun;
SSdHcan_all_lay_sun(isnan(SSdHcan_all_lay_sun)) = 0;

SSdHcan_all_lay = SSdHcan_all_lay_shade + SSdHcan_all_lay_sun;
SSdHcan_all = sum(SSdHcan_all_lay);        
SSdHcan = sum(SSdHcan_all);        


%**************************************************************************
%**************************************************************************
%                         II.  FOR THE SOIL
%**************************************************************************
%**************************************************************************
    % Latent Heat
    SSsoilLE_out = (LE_soil)./Tsurf;
    % Sensible Heat
    SSsoilH_out = (H_soil)./Tsurf;
    % Ground Heat Flux
    SSsoilG = G./Tsurf;
    %dH
    SSsoildH = dH./Tsurf;
    %dS
    SSsoildS = dS./Tsurf;
    
    
    
%**************************************************************************
%**************************************************************************
%**************************************************************************
%                         II.  THE NET ECOSYSTEM
%**************************************************************************
%**************************************************************************
%**************************************************************************

% canopy
SSLE_net =    (sum(sum(LE_can_all_lay))+LE_soil)/Tl_net;
LE_net =      (sum(sum(LE_can_all_lay))+LE_soil);
SSH_net =     (sum(sum(H_can_all_lay))+H_soil)/Tl_net;
H_net =       (sum(sum(H_can_all_lay))+H_soil);
dHcan_net =   (sum(sum(dHcan_all_lay)));
SSdHcan_net = (sum(sum(dHcan_all_lay)))/Tl_net;
SSG_net =     G/Tl_net;
G_net =       G;


%**************************************************************************
%**************************************************************************
%**************************************************************************
%               III.  STORE ENTROPY CALCULATIONS IN STRUCTURE
%**************************************************************************
%**************************************************************************
%**************************************************************************

% canopy
SSresults.SSLE_can_all_lay = SSLE_can_all_lay;                             %[W/m^2/K] Entropy due to LE in canopy for all layers and species
SSresults.SSLE_can_all = SSLE_can_all;                                     %[W/m^2/K] Entropy due to LE in canopy for all layers
SSresults.SSLE_can_tot = SSLE_can;                                         %[W/m^2/] Entropy due to LE 

SSresults.SSH_can_all_lay = SSH_can_all_lay;                               %[W/m^2/K] Entropy due to H in canopy for all layers and species 
SSresults.SSH_can_all = SSH_can_all;                                       %[W/m^2/K] Entropy due to H in canopy for all layers 
SSresults.SSH_can_tot = SSH_can;                                           %[W/m^2/] Entropy due to H 
SSresults.SSdHcan_tot = SSdHcan;                                           %[W/m^2/] Entropy due to change in heat in the leaves dHcan


% soil
SSresults.SSsoilLE_out = SSsoilLE_out;                                     %[W/m^2/K] Entropy due to LE in soil                                    
SSresults.SSsoilH_out = SSsoilH_out;                                       %[W/m^2/K] Entropy due to H in soil 
SSresults.SSsoilG = SSsoilG;                                               %[W/m^2/K] Entropy due to G in soil 
SSresults.SSsoildH = SSsoildH;                                             %[W/m^2/K] Entropy due to dH in soil 
SSresults.SSsoildS = SSsoildS;                                             %[W/m^2/K] Entropy due to dS in soil     

% net
SSresults.SSLE_net = SSLE_net;                                             %[W/m^2/K] Net (All Ecocystem) Entropy, due to LE 
SSresults.SSH_net = SSH_net;                                               %[W/m^2/K] Net (All Ecocystem) Entropy, due to H 
SSresults.SSG_net = SSG_net;                                               %[W/m^2/K] Net (All Ecocystem) Entropy, due to G 
SSresults.SSdHcan_net = SSdHcan_net;                                       %[W/m^2/K] Net (All Ecocystem) Entropy, due to dHcan 


% store net energy to compute Teffent
SSresults.LE_net = LE_net;                                                 %[W/m^2] Total (All Ecosystem) released LE
SSresults.H_net = H_net;                                                   %[W/m^2] Total (All Ecosystem) released H 
SSresults.dHcan_net = dHcan_net;                                           %[W/m^2] Total (All Ecosystem) dHcan  
SSresults.G_net = G_net;                                                   %[W/m^2] Total (All Ecosystem) G

