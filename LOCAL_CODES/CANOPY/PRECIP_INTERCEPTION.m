function [Sh2o_prof, Smaxz, ppt_ground, wetfrac, dryfrac] = PRECIP_INTERCEPTION(FORCING, VARIABLES, VERTSTRUC, PARAMS)

%=========================================================================
% Augments current canopy water storage with intercepted precipitation
%
% Written By: Darren Drewry
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING         % FORCING structure
%       VARIABLES       % VARIABLES structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%------------------------- Output Variables ------------------------------
%       Sh2o_prof       % [mm] Canopy moisture storage for all layers
%       Smaxz           % [mm] Maximum Canopy moisture storage 
%       ppt_ground      % [mm] Precipitation Reaching the ground
%       wetfrac         % [] Fraction of LAI that is wet
%       dryfrac         % [] Fraction of LAI that is dry
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % FORCING
    ppt = FORCING.ppt;                                                     % [mm] Precipitation Reaching the ground 
    
    % VARIABLES
    Sh2o_prof = VARIABLES.CANOPY.Sh2o_prof;                                % [mm] Canopy moisture storage for all layers 

    % VERTSTRUC
    LAIz = VERTSTRUC.LAIz;                                                 % [m^2 lea / m^2 ground]  LAI at all layers
    znc = VERTSTRUC.znc;                                                   % [] Height of canopy levels 
    nvinds = VERTSTRUC.nvinds;                                             % [] Indicator of Layers where LAI is zero of lower for all species  
    
    % PARAMS
    Smax = PARAMS.CanStruc.Smax;                                           % [mm] Maximum Canopy moisture storage  
    Ffact = PARAMS.CanStruc.Ffact;                                         % [] Max fraction of canopy that can be wet  
    pptintfact = PARAMS.CanStruc.pptintfact;                               % [] precip extinction coefficient
    clump = PARAMS.Rad.clump;                                              % [] Foliage clumping parameter 
%*************************************************************************
%*************************************************************************


    Smaxz = Smax * LAIz;    % [mm] Maximum canopy moisture storage
    
    % Attenuate precipitation through canopy
    if (ppt>0)

        % Loop from top to bottom of canopy
        for zz = length(znc):-1:1
            if (Sh2o_prof(zz) >= Smaxz(zz))
                continue;
            else
                % Potential H2O intercepted
                pptint_pot = (1 - exp(-pptintfact * clump * LAIz(zz))) * ppt;
                if ( (Sh2o_prof(zz) + pptint_pot) <= Smaxz(zz) )
                    Sh2o_prof(zz) = Sh2o_prof(zz) + pptint_pot;
                    ppt = ppt - pptint_pot;
                else
                    pptint = Smaxz(zz) - Sh2o_prof(zz);
                    Sh2o_prof(zz) = Smaxz(zz);
                    ppt = ppt - pptint;
                end
            end
            
        end
    end
    
    ppt_ground = ppt;  % ppt is now the precipitation incident on the ground
    
    wetfrac = Ffact.*(Sh2o_prof./Smaxz);
    % Adjust minimum fraction to 0.01
    wetfrac(wetfrac < 0.01) = 0;
    dryfrac = 1 - wetfrac;
    wetfrac(nvinds) = 0;
    dryfrac(nvinds) = 0;
