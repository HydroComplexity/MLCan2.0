function [fsun, fshade, LAIsun, LAIshade, ...
          SWabs_sun, SWabs_shade, PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, ...
          PARabs_sun_lai, PARabs_shade_lai, ...
          SWabs_soil, PARabs_soil, NIRabs_soil, ...
          SWout, PARout, NIRout, PARtop, NIRtop, fdiff, ...
          SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
          SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
          Kbm, taud, refl_soil] = ...
    SWRAD (FORCING, VERTSTRUC, PARAMS, CONSTANTS, VARIABLES)

%=========================================================================
% Augments current canopy water storage with intercepted precipitation
%
% Written By: Darren Drewry, Modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING          % VARIABLES structure
%       VERTSTRUC        % PARAMS structure
%       PARAMS           % CONSTANTS structure
%       CONSTANTS        % VARIABLES structure
%       VARIABLES        % VARIABLES structure
%------------------------- Output Variables ------------------------------
%       fsun            % [] Fraction of canopy that is sunlit
%       fshade          % [] Fraction of canopy that is shade
%       LAIsun          % [] Leaf Area Index that is sunlit
%       LAIshade        % [] Leaf Area Index that is shade
%       SWabs_sun       % [W/m^2 ground] Absorbed SW by sunlit fraction
%       SWabs_shade     % [W/m^2 ground] Absorbed SW by shade fraction
%       PARabs_sun      % [umol / m^2 ground / s] Absorbed PAR by sunlit fraction
%       PARabs_shade    % [umol / m^2 ground / s] Absorbed PAR shade by fraction
%       NIRabs_sun      % [W/m^2 ground] absorbed NIR by sunlit fraction 
%       NIRabs_shade    % [W/m^2 ground] absorbed NIR by shade fraction 
%       PARabs_sun_lai  % [umol / m^2 leaf / s] absorbed PAR per LAI by sunlit fraction  
%       PARabs_shade_lai% [umol / m^2 leaf / s] absorbed PAR per LAI by sunlit fraction  
%       SWabs_soil      % [W/m^2 ground] Absorbed SW by sunlit soil
%       PARabs_soil     % [umol / m^2 ground / s] Absorbed PAR by soil
%       NIRabs_soil     % [umol / m^2 ground / s] Absorbed NIR by soil 
%       SWout           % [W/m^2 ground] SW radiation out of the ecosystem
%       PARout          % [umol / m^2 ground / s] PAR out of the ecosystem
%       NIRout          % [W/m^2 ground] NIR out of the ecosystem
%       PARtop          % [umol / m^2 ground / s] Incident above-canopy PAR 
%       NIRtop          % [W/m^2 ground] Incident above-canopy NIR
%       fdiff           % [] Fraction of incident radiation that is diffused 
%       SWcandir_in     % [W/m^2 ground] Direct SW Radiation In, canopy 
%       SWcandir_out    % [W/m^2 ground] Direct SW Radiation Out, canopy
%       SWcandif_in     % [W/m^2 ground] Diffuse SW Radiation In, canopy
%       SWcandif_out    % [W/m^2 ground] Diffuse SW Radiation Out, canopy
%       SWsoildir_in    % [W/m^2 ground] Direct SW Radiation In, canopy
%       SWsoildir_out   % [W/m^2 ground] Direct SW Radiation Out, canopy 
%       SWsoildif_in    % [W/m^2 ground] Diffuse SW Radiation in canopy
%       SWsoildif_out   % [W/m^2 ground] Diffuse SW Radiation Out, canopy
%       Kbm             % [] Beam Extinction Coefficient 
%       taud            % [] Transmitted fraction
%       refl_soil       % [] Soil Reflectivity 
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************

% FORCING
SWin = FORCING.Rg;                                                         %[W/m^2] Incoming SW radiation   
zendeg = FORCING.zen;                                                      %[deg] zenith angle of sun  
doy = FORCING.doy;                                                         %[] Day of year 
Pa = FORCING.Pa;                                                           %[kPa] Air Pressure 

% VERTSTRUC
LAIz = VERTSTRUC.LAIz;                                                     %[m^2/m^2groun] LAI at every layer 
vinds = VERTSTRUC.vinds;                                                   %[] Indicator where LAI is lower or equal to zero  

% PARAMS
transmiss = PARAMS.Rad.transmiss;                                          %[] Atmospheric transmissivity  
xx = PARAMS.Rad.xx;                                                        %[] Leaf angle distribution param  
clump = PARAMS.Rad.clump;                                                  %[] Foliage clumping parameter 
trans_PAR = PARAMS.Rad.trans_PAR;                                          %[] Foliage transmissivity to PAR                                             
refl_PAR = PARAMS.Rad.refl_PAR;                                            %[] Foliage reflection PAR 
trans_NIR = PARAMS.Rad.trans_NIR;                                          %[] Foliage transmissivity to NIR 
refl_NIR = PARAMS.Rad.refl_NIR;                                            %[] Foliage reflection to NIR     
refl_snow = PARAMS.Rad.refl_snow;                                          %[] Reflectivity of new snow.  
refl_snow_old = PARAMS.Rad.refl_snow_old;                                  %[] How much reflectivity of snow decreases with age (per 12 days)
Kdf = PARAMS.Rad.Kdf;                                                      %[] Diffuse extinction coefficient  
    
% VARIABLES    
snow_tcount = VARIABLES.SOIL.snow_tcount;                                  %[] Number of time step with snow  
if VARIABLES.SOIL.zicesl > 5
    refl_soil = max((refl_snow - (snow_tcount/48/12.5)*refl_snow_old),0.1);%[] Soil Reflectivity
else
    %refl_soil = 0.31*exp(-12.7*volliq(1))+0.15;
    %refl_soil = 0.2;
    % Allison notices 2015.04.29
    refl_soil = PARAMS.Rad.refl_soil;
end
    
% CONSTANTS    
Wm2toumol = CONSTANTS.Wm2toumol;                                           %[] Conversion Factor from [W/m^2] to [umol]   
umoltoWm2 = CONSTANTS.umoltoWm2;                                           %[] Conversion Factor from [umol] to [W/m^2] 
fPAR = CONSTANTS.fPAR;                                                     %[] fraction incoming shortwave as PAR  
%*************************************************************************
%*************************************************************************


    zenrad = zendeg * pi/180;
%      if (zendeg > 89)
%          SWin = 0;
%          PARin = 0;
%      end
    
    PARtop = SWin*fPAR*Wm2toumol;
    NIRtop = SWin*(1-fPAR);

% Diffuse Fraction
    [fdiff] = DIFFUSE_FRACTION (zendeg, doy, SWin); 
    if (zendeg > 80)
        fdiff = 1;
    end
    
    PARtop_beam = PARtop * (1-fdiff);
    PARtop_diff = PARtop - PARtop_beam;
    NIRtop_beam = NIRtop * (1-fdiff);
    NIRtop_diff = NIRtop - NIRtop_beam;

    
% BEAM EXTINCTION COEFF FOR ELLIPSOIDAL LEAF DISTRIBUTION (EQN. 15.4 C&N)
    Kbm = sqrt(xx^2 + tan(zenrad)^2) / (xx + 1.774*(xx+1.182)^(-0.733));
        
% INITIALIZE ARRAYS
    len = length(LAIz);
    PARabs_sun = zeros(len,1);
    PARabs_shade = zeros(len,1);
    PARcandir_in = zeros(len,1);
    PARcandir_out = zeros(len,1);
    PARcandif_in = zeros(len,1);
    PARcandif_out = zeros(len,1);
    PARsoildir_in=0;
    PARsoildir_out=0;
    PARsoildif_in=0;
    PARsoildif_out=0;

    
    NIRabs_sun = zeros(len,1);
    NIRabs_shade = zeros(len,1);
    NIRcandir_in = zeros(len,1);
    NIRcandir_out = zeros(len,1);
    NIRcandif_in = zeros(len,1);
    NIRcandif_out = zeros(len,1);
    NIRsoildir_in=0;
    NIRsoildir_out=0;
    NIRsoildif_in=0;
    NIRsoildif_out=0;
      
    
    fsun = zeros(len,1);
    fshade = zeros(len,1);
    LAIsun = zeros(len,1);
    LAIshade = zeros(len,1);
    diffdn = zeros(len,1);
    diffup = zeros(len,1);
    PARabs_soil = 0;
    NIRabs_soil = 0;
    radlost = 0;
    PARout = 0;
    NIRout = 0;
    

if (PARtop<1)
    LAIshade = LAIz;
    fshade = fshade + 1;
else
    
    % Iterate to Solve PAR Absorption Profile
    count = 0;
    percdiff = 1;
    radin = PARtop; 
    beam_top = PARtop_beam; diff_top = PARtop_diff;
    while (percdiff > 0.01)

        [PARabs_sun, PARabs_shade, PARcandir_in, PARcandir_out, PARcandif_in, PARcandif_out,...
            PARabs_soil, PARsoildir_in, PARsoildir_out, PARsoildif_in, PARsoildif_out, fsun, fshade, diffdn, diffup, ...
                radabs_tot, radlost, PARradremain] = ...
            SW_ATTENUATION (beam_top, diff_top, LAIz, ...
                            trans_PAR, refl_PAR, refl_soil, clump, Kbm, Kdf, count, ...
                            PARabs_sun, PARabs_shade, PARcandir_in, PARcandir_out, PARcandif_in, PARcandif_out,...
                            diffdn, diffup, PARabs_soil, PARsoildir_in, PARsoildir_out, PARsoildif_in, PARsoildif_out, ...
                            fsun, fshade, radlost);
                        
         beam_top = 0;  diff_top = 0;

         radtot = radabs_tot + radlost;
         percdiff = (radin - radtot) / radin;

         count = count + 1;
         if (count>5)
             disp('COUNT > 5 in PAR loop!!!');
             break;
         end

    end
    
    LAIsun = fsun .* LAIz;
    LAIshade = fshade .* LAIz;
    
    PARout = radlost;


    % Iterate to Solve NIR Absorption Profile
    diffdn = zeros(len,1);
    diffup = zeros(len,1);
    radlost = 0;
    
    count = 0;
    percdiff = 1;
    radin = NIRtop; 
    beam_top = NIRtop_beam; diff_top = NIRtop_diff;
    while (percdiff > 0.01)

        [NIRabs_sun, NIRabs_shade, NIRcandir_in, NIRcandir_out, NIRcandif_in, NIRcandif_out,...
            NIRabs_soil, NIRsoildir_in, NIRsoildir_out, NIRsoildif_in, NIRsoildif_out, fsun, fshade, diffdn, diffup, ...
                radabs_tot, radlost, NIRradremain] = ...
            SW_ATTENUATION (beam_top, diff_top, LAIz, ...
                            trans_NIR, refl_NIR, refl_soil, clump, Kbm, Kdf, count, ...
                            NIRabs_sun, NIRabs_shade, NIRcandir_in, NIRcandir_out, NIRcandif_in, NIRcandif_out,...
                            diffdn, diffup, NIRabs_soil, NIRsoildir_in, NIRsoildir_out, NIRsoildif_in, NIRsoildif_out, ...
                            fsun, fshade, radlost);

         beam_top = 0;  diff_top = 0;

         radtot = radabs_tot + radlost;
         percdiff = (radin - radtot) / radin;

         count = count + 1;
         if (count>5)
             disp('COUNT > 5 in NIR loop!!!');
             break;
         end

    end
    NIRout = radlost;
        
end

taud = exp(-Kdf*clump*LAIz);

% Total SW outgoing [W/m^2]
    SWout = PARout*umoltoWm2 + NIRout;
    
% Shortwave Absorption Profiles    
    SWabs_sun = PARabs_sun*umoltoWm2 + NIRabs_sun;
    SWabs_shade = PARabs_shade*umoltoWm2 + NIRabs_shade;
    
% Shortwave Absorbed by Soil    
    SWabs_soil = NIRabs_soil + PARabs_soil*umoltoWm2;  % [W/m^2]
    
% Absorbed PAR per unit LAI [umol/m^2 LEAF AREA/s]
    PARabs_sun_lai = PARabs_sun;
    PARabs_sun_lai(vinds) = PARabs_sun_lai(vinds)./LAIsun(vinds);
    PARabs_sun_lai(LAIsun==0)=0;

    PARabs_shade_lai = PARabs_shade;
    PARabs_shade_lai(vinds) = PARabs_shade_lai(vinds)./LAIshade(vinds);
    PARabs_shade_lai(LAIsun==0) = 0 ;


% Shortwave Absorption. Direc and Diffuse
    SWcandir_in = PARcandir_in*umoltoWm2 + NIRcandir_in;
    SWcandir_out = PARcandir_out*umoltoWm2 + NIRcandir_out;
    SWcandif_in = PARcandif_in*umoltoWm2 + NIRcandif_in;
    SWcandif_out = PARcandif_out*umoltoWm2 + NIRcandif_out;    
       
    SWsoildir_in = PARsoildir_in*umoltoWm2 + NIRsoildir_in;
    SWsoildir_out = PARsoildir_out*umoltoWm2 + NIRsoildir_out;
    SWsoildif_in = PARsoildif_in*umoltoWm2 + NIRsoildif_in;
    SWsoildif_out = PARsoildif_out*umoltoWm2 + NIRsoildif_out;

    % Make output profiles column vectors
    fsun = fsun(:);
    fsun(fsun<1e-4) = 0;
    fshade = fshade(:);
    fshade(fshade<1e-4) = 0;    
    LAIsun = LAIsun(:);
    LAIsun(LAIsun<1e-4) = 0;
    LAIshade = LAIshade(:);
    LAIshade(LAIshade<1e-4) = 0;
    SWabs_sun = SWabs_sun(:);
    SWabs_sun(SWabs_sun<1e-4) = 0;
    SWabs_shade = SWabs_shade(:);
    SWabs_shade(SWabs_shade<1e-4) = 0;
    PARabs_sun = PARabs_sun(:);
    PARabs_sun(PARabs_sun<1e-4) = 0;
    PARabs_shade = PARabs_shade(:);
    PARabs_shade(PARabs_shade<1e-4) = 0;
    NIRabs_sun = NIRabs_sun(:);
    NIRabs_sun(NIRabs_sun<1e-4) = 0;
    NIRabs_shade = NIRabs_shade(:);
    NIRabs_shade(NIRabs_shade<1e-4) = 0;
  