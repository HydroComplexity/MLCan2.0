function [LWabs_can, LWabs_canM, LWabs_soil, LWabs_soilM, diffdn, diffdnM, diffup, diffupM, radlost, radlostM,...
    LW_can_out, LW_sun_out, LW_shade_out, LW_soil_out] = ...
        LW_ATTENUATION (LW_sky, Tatop, Tl_sun, Tl_shade, Tsoil, LAIz, ...
                       epsv, epss, clump, Kdf, count, fsun, fshade, ...
                       LAIsun, LAIshade, dz, LWabs_can, LWabs_canM, LWabs_soil, LWabs_soilM, ...
                       diffdn, diffdnM, diffup, diffupM, radlost, radlostM, LWmethod, retainLW)

%=========================================================================
%   Calls the longwave radiation absorption subroutine, iterating until most 
%       of the incident longwave, and that emitted by vegetation and the
%       soil, are absorbed or directed to the atmosphere
%
% Written By: Darren Drewry, Modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       LW_sky           % [W / m2 ground] Forcing LW radiation from atmosphere
%       Tatop            % [C] Temperature Atmosphere
%       Tl_sun           % [C] Temperature Canopy Sunlit Fraction
%       Tl_shade         % [C] Temperature Canopy Shade Fraction
%       Tsoil            % [C] Temperature Soil
%       LAIz             % [m^2 leaf / m^2 ground] LAI all canopy layers
%       epsv             % [] Vegetation emissivity
%       epss             % [] Soil emissivity
%       clump            % [] Foliage clumping parameter
%       kdf              % [] Diffuse extinction coefficient
%       count            % [] Counter longwave iteration
%       fsun             % [] Sunlit Fraction
%       fshade           % [] Shade Fraction
%       LAIsun           % [m^2 leaf / m^2 ground] LAI sunlit fraction
%       LAIshade         % [m^2 leaf / m^2 ground] LAI shade fraction
%       dz               % [m] Layer Thickness in the canopy
%       LWabs_can        % [m^2 leaf / m^2 ground] LW absorbed by canopy
%       LWabs_canM       % [W / m2 groun] Matrix of Absorptin of LW radiation in canopy, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWabs_soil       % [W / m^2 ground] absorbed LW by canopy in each layer
%       LWabs_soilM      % [W / m2 ground] Vector of Absorptin of LW radiation in the soil, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       diffdn           % [W / m2 ground] Downward radiation 
%       diffdnM          % [W / m2 ground] Matrix of downward radiation specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       diffup           % [W / m2 ground] Upward radiation 
%       diffupM          % [W / m2 ground] Matrix of upward radiation specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       radlost          % [W / m2 ground] Radiation Lost
%       radlostM         % [W / m2 ground] Vector of Radiation Lost specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWmethod         % [] Method to compute in LW_ATTENUATION FUNCITON LW 0. Original Darren. 1. Corrected
%       radlost          % [W / m2 ground] Radiation Lost
%       retainLW         % [] If LWmethod = 1,  retainLW is an extra LW retain the canopy
%------------------------- Output Variables ------------------------------
%       LWabs_can        % [W / m^2 ground] absorbed LW by canopy in each layer
%       LWabs_canM       % [W / m2 groun] Matrix of Absorptin of LW radiation in canopy, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWabs_soil       % [W / m^2 ground] absorbed LW by canopy in each layer
%       LWabs_soilM      % [W / m2 ground] Vector of Absorptin of LW radiation in the soil, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       diffdn           % [W / m2 ground] Downward radiation 
%       diffdnM          % [W / m2 ground] Matrix of downward radiation specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       diffup           % [W / m2 ground] Upward radiation 
%       diffupM          % [W / m2 ground] Matrix of upward radiation specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       radlost          % [W / m2 ground] Radiation Lost
%       radlostM         % [W / m2 ground] Vector of Radiation Lost specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LW_can_out       % [W / m2 ground] Radiation emitted by canopy 
%       LW_sun_out       % [W / m2 ground] Radiation emitted by canopy, sunlit fraction
%       LW_shade_out     % [W / m2 ground] Radiation emitted by canopy, shade fraction
%       LW_soil_out      % [W / m2 ground] Radiation emitted by soil
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
%
%   NOTE: index (1) refers to the layer just above the soil
%
%   Assume no reflected IR - all incident IR either absorbed or transmitted
%
%   Thermal emissivity and absorptivity are equal
%
%   Written By: Darren Drewry

soil_inc = 0;
soil_incM = 0;
boltz = 5.6697 * 10^-8;  % [W m^-2 K^-4]

%====================
% DOWNWARD RADIATION
%====================   
    tind = length(LAIz);
    bind = 1;
    diffdn(tind) = diffdn(tind) + LW_sky;
    diffdnM(tind,2*tind+1)= LW_sky;
    count2=0;
    for ii = tind:-1:bind 
        % ABSORBED DOWNWARD
            taud = exp(-Kdf * clump * LAIz(ii));

            LW_int = diffdn(ii) - taud*diffdn(ii);

            LWabs_can(ii) = LWabs_can(ii) + epsv*LW_int;
            
        % ****************************************************************                            
            changediff = (diffdnM(ii,:)./sum(diffdnM(ii,:)))*(epsv*LW_int);
            LWabs_canM(ii,:) = LWabs_canM(ii,:) + changediff;
        % ****************************************************************                            
        
        % DOWNWARD FLUX
            if (ii>bind)
                diffdn(ii-1) = diffdn(ii-1) + taud*diffdn(ii) + (1-epsv)*LW_int;
        % ****************************************************************                
                changediff = (diffdnM(ii,:)./sum(diffdnM(ii,:)))*(taud*diffdn(ii) + (1-epsv)*LW_int);
                diffdnM(ii-1,:) = diffdnM(ii-1,:) + changediff;
        % ****************************************************************                                
            else
                soil_inc = soil_inc + taud*diffdn(ii) + (1-epsv)*LW_int;
        % ****************************************************************                
                changediff = (diffdnM(ii,:)./sum(diffdnM(ii,:)))*(taud*diffdn(ii) + (1-epsv)*LW_int);
                soil_incM = soil_incM + changediff;
        % ****************************************************************                                          
            end
            
        % THERMAL CONTRIBUTION FROM FOLIAGE AND SOIL
            if (count==0)
                
                if LWmethod == 0
                     LW_flux_sun = fsun(ii)*(1-taud)*epsv*boltz*(Tl_sun(ii)+273.15)^4;
                     LW_flux_shade = fshade(ii)*(1-taud)*epsv*boltz*(Tl_shade(ii)+273.15)^4;
                     retainLW = 0;
                else
                     LW_flux_sun = (1-retainLW)*LAIsun(ii)*epsv*boltz*(Tl_sun(ii)+273.15)^4;
                     LW_flux_shade = (1-retainLW)*LAIshade(ii)*epsv*boltz*(Tl_shade(ii)+273.15)^4;                     
                end   
                                
                % Exclude nan values and set them to zero. nan means there
                % is not LAI at all in this layer, so no flux
                LW_flux_sun(isnan(LW_flux_sun))=0;
                LW_flux_shade(isnan(LW_flux_shade))=0;
                
                
                LW_flux = (LW_flux_sun + LW_flux_shade);                 
                delta = 2*(LW_flux/(1-retainLW) - LW_flux);
                
                LWabs_can(ii) = LWabs_can(ii) + delta;                                             
                LW_can_out(ii) = 2*LW_flux*1/(1-retainLW);
                LW_sun_out(ii) = 2*LW_flux_sun/(1-retainLW);
                LW_shade_out(ii) = 2*LW_flux_shade/(1-retainLW);
                
                
                
                % UPWARD FLUX
                if (ii<tind)
                    diffup(ii+1) = diffup(ii+1) + LW_flux;
         % ****************************************************************
                    diffupM(ii+1,ii) = diffupM(ii+1,ii) + LW_flux_sun;    % Emitted by sun part     
                    diffupM(ii+1,ii+tind) = diffupM(ii+1,ii+tind) + LW_flux_shade;    % Emitted by shade part
         % ****************************************************************                                                                                        
                else
                    radlost = radlost + LW_flux;
         % ****************************************************************
                    radlostM(1,ii) = LW_flux_sun;    % Emitted by sun part 
                    radlostM(1,ii+tind) = LW_flux_shade;    % Emitted by sun part 
         % ****************************************************************                
                end
                
                % DOWNWARD FLUX
                if (ii>bind)
                    diffdn(ii-1) = diffdn(ii-1) + LW_flux;
         % ****************************************************************                
                    diffdnM(ii-1,ii) = diffdnM(ii-1,ii) + LW_flux_sun;   % Emitted by sun part
                    diffdnM(ii-1,ii+tind) = diffdnM(ii-1,ii+tind) + LW_flux_shade;  % Emitted by shade part                  
         % ****************************************************************
                else
                   soil_inc = soil_inc + LW_flux;
         % ****************************************************************                
                   soil_incM(1,ii) = soil_incM(1,ii) + LW_flux_sun;     % Emitted by sun part 
                   soil_incM(1,ii+tind) = soil_incM(1,ii+tind) + LW_flux_shade;     % Emitted by shade part          
         % ****************************************************************                         
                end
                
                % SOIL EMISSION
                diffup(1) = epss*boltz*(Tsoil+273.15)^4; 
         % ****************************************************************                
                diffupM(1,2*tind+2) = diffup(1);
         % ****************************************************************                   
                LW_soil_out = epss*boltz*(Tsoil+273.15)^4;
            else
                LW_can_out = 0;
                LW_soil_out = 0;
            end
            
       diffdn(ii) = 0;  % Downward flux has been absorbed or transmitted
    end
    LWabs_soil = LWabs_soil + epss*soil_inc;

%   ***********************************************************************
    LWabs_soilM = LWabs_soilM + epss*soil_incM;
%   ***********************************************************************
    

    diffup(1) = diffup(1) + (1-epss)*soil_inc;
    % ****************************************************************                
    diffupM(1,:)=diffupM(1,:) + (soil_incM./sum(soil_incM))*(1-epss)*soil_inc;  
    % ****************************************************************
         
    
%==================
% UPWARD RADIATION
%==================      
    for ii = 1:tind
       
        % ABSORBED UPWARD
            taud = exp(-Kdf * clump * LAIz(ii));

            LW_int = diffup(ii) - taud*diffup(ii);

            LWabs_can(ii) = LWabs_can(ii) + epsv*LW_int;
            
        % ****************************************************************                            
            changediff = (diffupM(ii,:)./sum(diffupM(ii,:)))*(epsv*LW_int);
            LWabs_canM(ii,:) = LWabs_canM(ii,:) + changediff;
        % ****************************************************************                            

            
        % UPWARD FLUX
            if (ii<tind)
                diffup(ii+1) = diffup(ii+1) + taud*diffup(ii) + (1-epsv)*LW_int;
         % ****************************************************************                
                changediff = (diffupM(ii,:)./sum(diffupM(ii,:)))*(taud*diffup(ii) + (1-epsv)*LW_int);
                diffupM(ii+1,:) = diffupM(ii+1,:) + changediff;
         % ****************************************************************                    
            else
                radlost = radlost + taud*diffup(ii) + (1-epsv)*LW_int;
         % ****************************************************************                
                changediff = (diffupM(ii,:)./sum(diffupM(ii,:)))*(taud*diffup(ii) + (1-epsv)*LW_int);
                radlostM = radlostM + changediff;
         % ****************************************************************                
            end
       
        diffup(ii) = 0;
    end
    
