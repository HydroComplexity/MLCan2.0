function [LWabs_can, LWabs_canM, LWabs_sun, LWabs_shade, LWabs_soil, LWabs_soilM, LW_sky, LW_sky2, LWup, LWupM, ...
          LWcan_emit, LWsun_emit, LWshade_emit, LWsoil_emit] = ...
    LWRAD (FORCING, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS)

%=========================================================================
%   Calls the longwave radiation absorption subroutine, iterating until most 
%       of the incident longwave, and that emitted by vegetation and the
%       soil, are absorbed or directed to the atmosphere
%
% Written By: Darren Drewry, Modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING          % VARIABLES structure
%       VARIABLES        % VARIABLES structure
%       VERTSTRUC        % PARAMS structure
%       PARAMS           % CONSTANTS structure
%       CONSTANTS        % VARIABLES structure
%------------------------- Output Variables ------------------------------
%       LWabs_can        % [W / m^2 ground] absorbed LW by canopy in each layer
%       LWabs_canM       % [W / m2 groun] Matrix of Absorptin of LW radiation in canopy, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWabs_sun        % [W / m^2 ground] absorbed LW by canopy sunlit fraction in each layer
%       LWabs_shade      % [W / m^2 ground] absorbed LW by canopy sunlit fraction in each layer
%       LWabs_soil       % [W / m^2 ground] absorbed LW by canopy in each layer
%       LWabs_soilM      % [W / m2 ground] Vector of Absorptin of LW radiation in the soil, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LW_sky           % [W / m^2 ground] Downward LW radiation 1
%       LW_sky2          % [W / m^2 ground] Downward LW radiation 2
%       LWup             % [W / m^2 ground] total leaving LW from canopy
%       LWupM            % [W / m^2 ground] Vector of leaving LW from with information from LW emmited from different layers and atmosphere
%       LWcan_emit       % [W / m^2 ground] LW emitted by canopy
%       LWcan_emit       % [W / m^2 ground] LW emitted by canopy sunlit fraction
%       LWcan_emit       % [W / m^2 ground] LW emitted by canopy shade fraction
%       LWsoil_emit      % [W / m^2 ground] LW emitted from soil
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************

% FORCING
zendeg = FORCING.zen;                                                      % [deg] zenith angle of sun  
Tatop = FORCING.Ta;                                                        % [C] Atmospheric Temperature 
eatop = FORCING.ea;                                                        % [kPa] Atmospheric Vapor Pressure
%    ELEV=FORCING.ELEV;
%    HR=FORCING.HR;
%    Pa=FORCING.Pa;
    
% VARIABLES
Tl_sun = VARIABLES.CANOPY.Tl_can_sun;                                      % [C] Canopy Temperature, Sunlit Fraction 
Tl_shade = VARIABLES.CANOPY.Tl_can_shade;                                  % [C] Canopy Temperature, Shade Fraction 
fsun = VARIABLES.CANOPY.fsun;                                              % [] Sunlit Fraction 
fshade = VARIABLES.CANOPY.fshade;                                          % [] Shade Fraction
LAIsun = VARIABLES.CANOPY.LAIsun;                                          % [m^2 leaaf/ m^2 ground] LAI sunlit             
LAIshade = VARIABLES.CANOPY.LAIshade;                                      % [m^2 leaaf/ m^2 ground] LAI shade     
Ts = VARIABLES.SOIL.Tsurf;
 
% VERTSTRUC
LAIz = VERTSTRUC.LAIz;                                                     % [m^2 leaaf/ m^2 ground] LAI at every layer 
dz = VERTSTRUC.dzc;                                                        % [m] Layer Thickness in the canopy 
    
% PARAMS
Kdf = PARAMS.Rad.Kdf;                                                      % [] Diffuse extinction coefficient
clump = PARAMS.Rad.clump;                                                  % [] Foliage clumping parameter 
epsv = PARAMS.Rad.epsv;                                                    % [] Vegetation emissivity 
epss = PARAMS.Rad.epss;                                                    % [] Soil emissivity 
LWcom = PARAMS.LWcom;                                                      % [] Computation of longwave incoming 1.DATA 2.Boltzman 3.Boltzman Corrected  
LWmethod = PARAMS.LWmethod;                                                % [] Method to compute in LW_ATTENUATION FUNCITON LW 0. Original Darren. 1. Corrected 
retainLW = PARAMS.retainLW;                                                % [] If LWmethod = 1,  retainLW is an extra LW retain the canopy

% CONSTANTS
boltz = CONSTANTS.boltz;                                                   % [W/m^2/K^4] Stefan-Boltzmann constant 
 
%*************************************************************************
%*************************************************************************


zenrad = zendeg * pi/180;
    
% DOWNWARD LW FROM SKY
   if LWcom ==1  % Using Longwave Incoming Data
        LWin = FORCING.LWdn;
        LW_sky = LWin;
        LW_sky2 = nan;  % use LW_sky2 in order to compare different methods
        if isnan(LW_sky)
            epsa = 1.72 * (eatop / (Tatop+273.15))^(1/7);
            LW_sky = epsa*boltz*(Tatop+273.15)^4;
        end                    
   elseif LWcom == 2; % Using Boltzman  
        epsa = 1.72 * (eatop / (Tatop+273.15))^(1/7);
        LW_sky = epsa*boltz*(Tatop+273.15)^4;
        LW_sky2 = nan; 
   elseif  LWcom == 3  % Using Boltzman corrected  
        epsa = 1.72 * (eatop / (Tatop+273.15))^(1/7);
        LW_sky2 = epsa*boltz*(Tatop+273.15)^4; 
        c=0.5;
        epsac=epsa;
%        epsac=9.2*10^(-6)*(Tatop+273.15)^(2);
        epsa=(1-0.84*c)*epsac+0.84*c;
        LW_sky = epsa*boltz*(Tatop+273.15)^4;
   end

%    %********* Computing longwave with Marks and Dozier (1979)     
%         Tcor=(Tatop+273.15)+(0.0065*ELEV);        
%         aa = 0.611;  % [kPa]
%         bb = 17.502;
%         cc = 240.97; % [C]
%         esatcor = aa*exp((bb*(Tcor-273.15))./(cc+(Tcor-273.15)));
%         ecor=HR*esatcor;    
%         epsa = (1.72 * (ecor/Tcor)^(1/7))*(Pa/101.3);
%         LW_sky2 = 0.9*epsa*boltz*(Tatop+273.15)^4+0.1*epss*boltz*(Ts+273.15)^4;        
   %********* Computing longwave with Idso and Jackson (1969)
%         LW_sky=(1-0.261*exp(-7.77*10^(-4)*(Tatop)^(2)))*boltz*(Tatop+273.15)^4;
   %********* Computing longwave with Idso (1974)
   %      LW_sky2=(1-0.261*exp(-7.77*10^(-4)*(Tatop)^(2)))*boltz*(Tatop+273.15)^4;
         
     
%     end

    
% INITIALIZE ARRAYS
    len = length(LAIz);
    LWabs_can = zeros(len,1);
    LWabs_canM = zeros(len,2*len+2);
    diffdn = zeros(len,1);
    diffdnM = zeros(len,2*len+2);
    diffup = zeros(len,1);
    diffupM = zeros(len,2*len+2);
    LWabs_soil = 0;
    LWabs_soilM = 0;
    radlost = 0;
    radlostM = zeros(1,2*len+2);
    
% Iterate to Solve LW Absorption Profile
    count = 0;
    maxiters = 10;
    percdiff = 1;
    LW_top = LW_sky; 
                     
    while (percdiff > 0.01)
        
        [LWabs_can, LWabs_canM, LWabs_soil, LWabs_soilM, diffdn, diffdnM, diffup, diffupM, radlost, radlostM, LW_can_out, LW_sun_out, LW_shade_out, LW_soil_out] = ...
            LW_ATTENUATION (LW_top, Tatop, Tl_sun, Tl_shade, Ts, LAIz, ...
                           epsv, epss, clump, Kdf, count, fsun, fshade, ...
                           LAIsun, LAIshade, dz, LWabs_can, LWabs_canM, LWabs_soil, LWabs_soilM,...
                           diffdn,diffdnM, diffup, diffupM, radlost, radlostM, LWmethod, retainLW);
                            
        if (count==0)
           LWtot = LW_sky + sum(LW_can_out) + LW_soil_out; 
           
           LWcan_emit = LW_can_out(:);
           
           LWsun_emit = LW_sun_out(:);
           LWshade_emit = LW_shade_out(:);
           
           LWsoil_emit = LW_soil_out;
           
        end

        LW_top = 0;  

        radremain = sum(diffdn) + sum(diffup);
        percdiff = radremain / LWtot;

        count = count + 1;
        if (count>maxiters)
            disp('TOO MANY ITERATIONS IN  LW LOOP!!!');
            break;
        end

    end
   
    LWup = radlost;
    LWupM = radlostM;
    
% Total Absorbed Radiation By Each Canopy Fraction (Sunlit and Shaded)
    LWabs_sun = LWabs_can.*fsun;
    LWabs_shade = LWabs_can.*fshade;