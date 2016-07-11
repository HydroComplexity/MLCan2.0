function [Sh2o_prof, Sh2o_can, VARIABLES, Evap_can, Evap_prof, ppt_ground] = EVAP_CONDENSATION_ADJUST(VARIABLES, VERTSTRUC, PARAMS)

%=========================================================================
% Adjust canopy water storage for condensation and evaporation
%
% Written By: Darren Drewry
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES       % VARIABLES structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%------------------------- Output Variables ------------------------------
%       Sh2o_prof       % [mm] Canopy moisture storage for all layers
%       Sh2o_can        % [mm] Canopy moisture storage for all layers
%       VARIABLES       % VARIABLES structure
%       Evap_can        % [mm] Evaporation from canopy 
%       Evap_prof       % [mm] Evaporation from canopy at all layers
%       ppt_ground      % [mm] Precipitation Reaching the ground
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************

%   VARIABLES
    Ch2o_prof = VARIABLES.CANOPY.Ch2o_prof;                                % [mm / s / m^2 leaf] Condensation water in the canopy 
    Evap_prof = VARIABLES.CANOPY.Evap_prof;                                % [mm] Evaporation from canopy at all layers  
    Sh2o_prof = VARIABLES.CANOPY.Sh2o_prof;                                % [mm] Canopy moisture storage for all layers 
    Sh2o_prof_ini = Sh2o_prof;                                             % [mm] Canopy moisture storage for all layers initial 
    Smaxz = VARIABLES.CANOPY.Smaxz;                                        % [mm] Maximum canopy moisture storage     
    ppt_ground = VARIABLES.SOIL.ppt_ground;                                % [mm] Precipitation Reaching the ground
    
%   VERTSTRUC    
    znc = VERTSTRUC.znc;                                                   % [] Height of canopy levels     

%   PARAMS    
    Ffact = PARAMS.CanStruc.Ffact;                                         % [] Max fraction of canopy that can be wet 
%*************************************************************************
%*************************************************************************
    

% Adjust Canopy Water Storage
        H2oinc = Ch2o_prof - Evap_prof;
        drip(length(znc)) = 0;
        % Loop from top to bottom of canopy
        dripout = 0;
        for zz = length(znc):-1:1
            Sh2o_prof(zz) = Sh2o_prof(zz) + H2oinc(zz) + drip(zz);            
            if (Sh2o_prof(zz) < 0)
                Evap_prof(zz) = Evap_prof(zz) + Sh2o_prof(zz); 
                Sh2o_prof(zz) = 0;
            elseif (Sh2o_prof(zz) >= Smaxz(zz))
                if zz==1
                    dripout =  Sh2o_prof(zz) - Smaxz(zz);
                else
                    drip(zz-1) = Sh2o_prof(zz) - Smaxz(zz);
                end
                Sh2o_prof(zz) = Smaxz(zz);
            end
        end
        Sh2o_can = sum(Sh2o_prof);
        ppt_ground = ppt_ground + dripout;  % total H2O incident on the ground
        wetfrac = Ffact.*(Sh2o_prof./Smaxz);
        dryfrac = 1 - wetfrac;
        
        % Modify Evaporation (Can no be higher than available storage)    
        VARIABLES.CANOPY.Evap_prof = Evap_prof;        
        Evap_can = sum(Evap_prof); % [mm] Total canopy evaporation

        %ASSIGN
        VARIABLES.SOIL.ppt_ground = ppt_ground;
        VARIABLES.CANOPY.Evap_prof = Evap_prof;

        
 %   error = sum(Ch2o_prof) - Evap_can - sum(Sh2o_prof - Sh2o_prof_ini) - ppt_ground; 
        
        
        
        
        
        
        
        