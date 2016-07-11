

function [sun_abs, shade_abs, candir_in, candir_out, candif_in, candif_out,...
    soil_abs, soildir_in, soildir_out, soildif_in, soildif_out, fsun, fshade, diffdn, diffup, ...
            radabs_tot, radlost, radremain] = ...
        SW_ATTENUATION (beam_top, diff_top, LAIz, ...
                        trans, refl, refl_soil, clump, Kbm, Kdf, count, ...
                        sun_abs, shade_abs, candir_in, candir_out, candif_in, candif_out ,...
                        diffdn, diffup, soil_abs, soildir_in, soildir_out, soildif_in, soildif_out, fsun, fshade, ...
                        radlost)
                    
%=========================================================================
% Augments current canopy water storage with intercepted precipitation
%
% Written By: Darren Drewry, modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       beam_top         % [umol / m^2 ground / s] Beam PAR at the top 
%       diff_top         % [umol / m^2 ground / s] Diffuse PAR at the top 
%       LAIz             % [m^2 leaf/m^2 ground] Leaf Area Index at each layer
%       trans            % [] Foliage transmissivity to PAR
%       refl             % [] Reflection to PAR canopy
%       refl_soil        % [] Soil Reflection
%       clump            % [] Foliage clumping parameter 
%       Kbm              % [] Beam Extinction Coefficient
%       Kdf              % [] Diffuse extinction coefficient
%       count            % [] Number of time step with snow
%       sun_abs          % [W/m^2 ground] Absorbed SW by sunlit fraction
%       shade_abs        % [W/m^2 ground] Absorbed SW by shade fraction
%       candir_in        % [W/m^2 ground] Direct SW Radiation In, canopy
%       candir_out       % [W/m^2 ground] Direct SW Radiation Out, canopy
%       candif_in        % [W/m^2 ground] Diffuse SW Radiation In, canopy
%       candif_out       % [W/m^2 ground] Diffuse SW Radiation Out, canopy
%       diffdn           % [W/m^2 ground] Diffuse SW Radiation down
%       diffup           % [W/m^2 ground] Diffuse SW Radiation up
%       soil_abs         % [W/m^2 ground] Soil absorption
%       soildir_in       % [W/m^2 ground] Direct SW Radiation In, soil
%       soildir_out      % [W/m^2 ground] Direct SW Radiation In, soil
%       soildif_in       % [W/m^2 ground] Direct SW Radiation In, soil
%       soildif_out      % [W/m^2 ground] Direct SW Radiation In, soil
%       fsun             % [] Radiation fraction sunlit
%       fshade           % [] Radiation fraction shade
%       radlost          % [W/m^2 ground] Radiation Lost

%------------------------- Output Variables ------------------------------
%       sun_abs          % [W/m^2 ground] Absorbed SW by sunlit fraction
%       shade_abs        % [W/m^2 ground] Absorbed SW by shade fraction
%       candir_in        % [W/m^2 ground] Direct SW Radiation In, canopy
%       candir_out       % [W/m^2 ground] Direct SW Radiation Out, canopy
%       candif_in        % [W/m^2 ground] Diffuse SW Radiation In, canopy
%       candif_out       % [W/m^2 ground] Diffuse SW Radiation Out, canopy
%       soil_abs         % [W/m^2 ground] Soil absorption
%       soildir_in       % [W/m^2 ground] Direct SW Radiation In, soil
%       soildir_out      % [W/m^2 ground] Direct SW Radiation In, soil
%       soildif_in       % [W/m^2 ground] Direct SW Radiation In, soil
%       soildif_out      % [W/m^2 ground] Direct SW Radiation In, soil
%       fsun             % [] Radiation fraction sunlit
%       fshade           % [] Radiation fraction shade
%       diffdn           % [W/m^2 ground] Diffuse SW Radiation down
%       diffup           % [W/m^2 ground] Diffuse SW Radiation up
%       radabs_tot       % [W/m^2 ground] Total absrbed radiation
%       radlost          % [W/m^2 ground] Radiation Lost
%       radremain        % Radiation Remaining in the System
% 
%========================================================================
                    
                    
%
%   NOTE: index (1) refers to the layer just above the soil
%

soil_inc = 0;

%====================
% DOWNWARD RADIATION
%====================   
    tind = length(LAIz);
    bind = 1;
    beam_inc = beam_top;
    diffdn(tind) = diffdn(tind) + diff_top;
    LAIc = 0;
    
    soildif_inc = 0;
    soildir_inc = 0;
    for ii = tind:-1:bind
        
        LAIc = LAIc + LAIz(ii);
        fsun(ii) = exp(-Kbm * clump * LAIc); % sunlit leaf fraction
        fshade(ii) = 1 - fsun(ii);    
            
      % BEAM RADIATION     
        if (Kbm>0 && beam_top>0 && count==0)
   
            taub = exp(-Kbm * clump * LAIz(ii));

            beam_int = beam_inc - taub*beam_inc; % intercepted beam

            beam_inc = taub*beam_inc;            % incident radiation on layer below

            sun_abs(ii) = (1-refl-trans)*beam_int;
            % ******************* to compute entropy *********************    
            candir_in(ii) = candir_in(ii) + beam_int;
            candir_out(ii) = candir_out(ii) + beam_int*(trans);
            candif_in(ii) = candif_in(ii) + 0;
            candif_out(ii) = candif_out(ii) +  beam_int*(refl);            
            % ******************* to compute entropy *********************    
            % intercepted beam that is transmitted
            if (ii==bind)
                soil_inc = trans*beam_int + beam_inc;
                soildir_inc = soildir_inc + trans*beam_int + beam_inc;
            else
                diffdn(ii-1) = trans*beam_int;
            end
            
            % intercepted beam that is reflected
            if (ii<tind)
                diffup(ii+1) = refl*beam_int;
            else
                reflbup = refl*beam_int; % used for debugging
            end
        else
            reflbup = 0;
        end
 
        
      % DIFFUSE RADIATION
        taud = exp(-Kdf * clump * LAIz(ii));

        diff_int = diffdn(ii) - taud*diffdn(ii); % intercepted downward diffuse
        
        % downward transmission
        if (ii==bind)
            soil_inc = soil_inc + trans*diff_int + taud*diffdn(ii);
            soildif_inc = soildif_inc + trans*diff_int + taud*diffdn(ii);
        else
            diffdn(ii-1) = diffdn(ii-1) + trans*diff_int + taud*diffdn(ii);
        end
        
        % upward reflection
        if (ii<tind)
            diffup(ii+1) = diffup(ii+1) + refl*diff_int;
        else
            refldup = refl*diff_int; % for debugging
        end

        % absorbed fraction
        sun_abs(ii) = sun_abs(ii) + (1-refl-trans)*diff_int*fsun(ii);
        shade_abs(ii) = shade_abs(ii) + (1-refl-trans)*diff_int*fshade(ii);
        
        % ******************* to compute entropy *********************    
        candif_in(ii) = candif_in(ii) + diff_int;
        candif_out(ii) = candif_out(ii) + diff_int*(refl+trans); 
        % ******************* to compute entropy *********************    
       
    end
    soil_abs = soil_abs + (1-refl_soil)*soil_inc;
    % ******************* to compute entropy *********************    
    soildir_in = soildir_in + soildir_inc;
    soildir_out = soildir_out + 0;
    soildif_in = soildif_in + soildif_inc;
    soildif_out = soildif_out + soildir_inc*refl_soil + soildif_inc*refl_soil;
    % ******************* to compute entropy *********************    
                
    diffup(1) = diffup(1) + refl_soil*soil_inc;%diffdn(1);
    soil_inc = 0;
    soildif_inc = 0;
    soildir_inc = 0;
    
    diffdn = diffdn*0;
    
    
%==========================
% UPWARD DIFFUSE RADIATION
%==========================
    for ii = bind:tind
        taud = exp(-Kdf * clump * LAIz(ii));
        
        diff_int = diffup(ii) - taud*diffup(ii);
        
        % upward tranmission
        if (ii<tind)
            diffup(ii+1) = diffup(ii+1) + trans*diff_int + taud*diffup(ii);
        else
            diffuplost = trans*diff_int + taud*diffup(ii); 
        end
        
        % downward reflection
        if (ii==bind)
            soil_inc = soil_inc + refl*diff_int;
            soildif_inc = soildif_inc + refl*diff_int;            
        else
            diffdn(ii-1) = diffdn(ii-1) + refl*diff_int;
        end
        
        % absorbed fraction
        sun_abs(ii) = sun_abs(ii) + (1-refl-trans)*diff_int*fsun(ii);
        shade_abs(ii) = shade_abs(ii) + (1-refl-trans)*diff_int*fshade(ii);
    % ******************* to compute entropy *********************    
        candif_in(ii) = candif_in(ii) + diff_int;
        candif_out(ii) = candif_out(ii) + diff_int*(refl+trans);
    % ******************* to compute entropy *********************    
    end
    diffup = diffup*0;
    soil_abs = soil_abs + (1-refl_soil)*soil_inc;

    soildif_in = soildif_in + soildif_inc;
    soildif_out = soildif_out + soildif_inc*(refl_soil);
            
    diffup(1) = refl_soil*soil_inc;
    
    
% Absorbed Radiation
    radabs_tot = sum(sun_abs) + sum(shade_abs) + soil_abs;
% Lost Radiation
    radlost = radlost + reflbup + refldup + diffuplost;
% Radiation Remaining in System
    radremain = sum(diffdn) + sum(diffup);