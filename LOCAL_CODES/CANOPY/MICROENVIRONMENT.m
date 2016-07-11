function [CAz, EAz, TAz] = MICROENVIRONMENT(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS)

%=========================================================================
% Calculate the canopy microenvironment (Ca, Ta, ea) using a first-order
% canopy closure model
%
% Written By: Darren Drewry
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING         % FORCING structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%       CONSTANTS       % CONSTANTS structure
%------------------------- Output Variables ------------------------------
%       CAz             %[umol/mol] Atmospheric Concentration of CO2, for all layers
%       EAz             %[kPa] Vapor Pressure Air, for all layers
%       TAz             %[C] Air Temperature, for all layers
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
% Calculate the canopy microenvironment (Ca, Ta, ea) using a first-order
% canopy closure model
%
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES
    CAz = VARIABLES.CANOPY.CAz;                                            % [umol/mol] Atmospheric Concentration of CO2, for all layers 
    TAz = VARIABLES.CANOPY.TAz;                                            % [C] Air Temperature, for all layers 
    EAz = VARIABLES.CANOPY.EAz;                                            % [kPa] Vapor Pressure Air, for all layers 
    PAz = VARIABLES.CANOPY.PAz;                                            % [kPa] Pressure Air, for all layers  
    Km = VARIABLES.CANOPY.Km;                                              % Momentum Diffusivity Constant      
    An_sun = VARIABLES.CANOPY.An_sun;                                      % [umol CO2/ m^2 leaf / s] Photosynthetic minus Leaf Respiration flux from canopy sunlit  
    An_shade = VARIABLES.CANOPY.An_shade;                                  % [umol CO2/ m^2 leaf / s] Photosynthetic minus Leaf Respiration flux from canopy shade
    LE_sun = VARIABLES.CANOPY.LE_sun;                                      % [W / m^2 leaf] LE flux from canopy sunlit     
    LE_shade = VARIABLES.CANOPY.LE_shade;                                  % [W / m^2 leaf] LE flux from canopy shade  
    H_sun = VARIABLES.CANOPY.H_sun;                                        % [W / m^2 leaf] H flux from canopy shade     
    H_shade = VARIABLES.CANOPY.H_shade;                                    % [W / m^2 leaf] H flux from canopy shade 
    LAIsun = VARIABLES.CANOPY.LAIsun;                                      % [m^2 leaf / m^2 ground] LAI sunlit
    LAIshade = VARIABLES.CANOPY.LAIshade;                                  % [m^2 leaf / m^2 ground] LAI shade
    Fc_soil = VARIABLES.SOIL.Fc_soil;                                      % [m^2 leaf / m^2 ground] CO2 flux from soil
    LE_soil = VARIABLES.SOIL.LE_soil;                                      % [m^2 leaf / m^2 ground] LE flux from soil 
    H_soil = VARIABLES.SOIL.H_soil;                                        % [m^2 leaf / m^2 ground] H flux from soil 

    % VERTSTRUC
    znc = VERTSTRUC.znc;                                                   % [m] Height of canopy levels 
    dzc = VERTSTRUC.dzc;                                                   % [m] Tickness of layers in canopy  
    fLAIz = VERTSTRUC.fLAIz;                                               % [] Fraction of total LAI for each species at each canopy layer    
            
    % PARAMS
    hcan = PARAMS.CanStruc.hcan;                                           % [h] Canopy Height  
    
    % CONSTANTS
    Lv = CONSTANTS.Lv;                                                     % Latent heat of vaporization 
    cp_mol = CONSTANTS.cp_mol;                                             %[J/mol/K] specific heat of air at constant pressure      
    
%*************************************************************************
%*************************************************************************

    molar_density = 44.6 * PAz * 273.15 ./ (101.3 * (TAz + 273.15));
    psy = 6.66 * 10^-4; % [1/C]

     % CO2
        %Ca = Ca .* molar_density;  % [umol / m^3]
        Sc = ((sum(An_sun.*fLAIz,2).*LAIsun) + (nansum(An_shade.*fLAIz,2).*LAIshade))./dzc;
        Sc = -Sc./molar_density;        % [umol/mol / s]
        [CAz] = ORDER_1_CLOSURE_ALL ( CAz, znc, dzc, Km, Sc, Fc_soil, hcan );  

     % VAPOR
        q = (EAz./PAz).*molar_density;                 
        Sv = ((sum(LE_sun.*fLAIz,2).*LAIsun) + (sum(LE_shade.*fLAIz,2).*LAIshade))./dzc; % [W / m^3]
        Sv = (Sv./Lv);
        Sv_soil =  LE_soil/Lv;
        [q] = ORDER_1_CLOSURE_ALL ( q, znc, dzc, Km, Sv, Sv_soil, hcan );
        EAz = (q ./ molar_density) .* PAz;

     % HEAT
        heat = TAz .* molar_density .* psy;
        Sh = ((sum(H_sun.*fLAIz,2).*LAIsun) + (sum(H_shade.*fLAIz,2).*LAIshade))./dzc; % [W / m^3]        
        Sh = Sh./cp_mol./molar_density;
        Sh_soil = H_soil./cp_mol./molar_density(1);
        [TAz] = ORDER_1_CLOSURE_ALL ( TAz, znc, dzc, Km, Sh, Sh_soil, hcan );   
        
        