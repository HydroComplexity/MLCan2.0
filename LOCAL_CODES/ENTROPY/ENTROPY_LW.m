
function [SSresults] = ENTROPY_LW (LWabs_canM, LWabs_soilM, LWemit_soil, LWemit_sun, LWemit_shade,...
                                   LWin, LWout, Tatop, Tsurf, boltz, entropymethod, zicesl, ...
                                   PARAMS, VARIABLES, CONSTANTS, SWITCHES, SSresults)

%=========================================================================
% This code computes the fluxes of entropy due to longwave. 
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       LWabs_canM                     % [W/m2] Matrix of Absorptin of LW radiation in canopy, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWabs_soilM                    % [W/m2] Vector of Absorptin of LW radiation in the soil, specifying the source from all the canopy layers (sunlit and shade), soil, and atmosphere
%       LWemit_soil                    % [W/m2] LW radiation emitted from soil
%       LWemit_sun                     % [W/m2] LW radiation emitted from canopy sunlit
%       LWemit_shade                   % [W/m2] LW radiation emitted from canopy sunlit
%       LWin                           % [W/m2] Total LW radiation in 
%       LWout                          % [W/m2] Total LW radiation out 
%       Tatop                          % [C] Atmospheric Temperature
%       Tsurf                          % [C] Soil Surface Temperature
%       boltz                          % [W/m^2/K^4] Stefan-Boltzmann constant
%       entropymethod                  % [] Method used to compute the entropy in radiation
%       zicesl                         % [mm] Ice water depth in Snow-litter pack 
%       PARAMS                         % PARAMS structure
%       VARIABLES                      % VARIABLES structure
%       CONSTANTS                      % CONSTANTS structure
%       SWITCHES                       % SWITCHES structure
%       SSresults                      % SSresults structure 
%------------------------- Output Variables ------------------------------
%       SSresults                      % SSresults structure 
% 
%========================================================================
% De reference block
%*********************************************************************
             
% PARAMS
epsv=PARAMS.Rad.epsv;                                                      %[] vegetation emissivity 
epss=PARAMS.Rad.epss;                                                      %[] soil emissivity 
epsa=PARAMS.Rad.epsa;                                                      %[] atmospheric emissivity
c2 = PARAMS.Entropy.c2;                                                    %[] Constant to conpute entropy radiation 
c3 = PARAMS.Entropy.c3;                                                    %[] Constant to conpute entropy radiation 
LWmethod = PARAMS.LWmethod;                                                %[] Method used to compute LW in canopy 
if LWmethod == 0
   retainLW = 0;
else 
   retainLW = PARAMS.retainLW;
end

% VARIABLES
LAIsun = VARIABLES.CANOPY.LAIsun;                                          %[m^2 leaf / m^2 ground] LAI in sunlit 
LAIshade = VARIABLES.CANOPY.LAIshade;                                      %[m^2 leaf / m^2 ground] LAI in shade 
fsun = VARIABLES.CANOPY.fsun;                                              %[] fraction sunlit 
fshade = VARIABLES.CANOPY.fshade;                                          %[] fraction shade 
plants = SWITCHES.plants;                                                  %[] indicator, plants or no plants at all     

% CONSTANTS
bol = CONSTANTS.boltz;

LAI = nansum(LAIsun + LAIshade);
%**************************************************************************
% compute number of layers
nlayers=PARAMS.CanStruc.nl_can; 

%**************************************************************************
% ALLOCATE MEMORY
SScanLW_in = zeros(nlayers,1);
%SScanLW_in = zeros(nlayers,1);

%**************************************************************************
% change temperature units from C to K
Tatop = Tatop + 273.15;
Tsurf = Tsurf + 273.15;

Tl_can_sun = VARIABLES.CANOPY.Tl_can_sun + 273;
Tl_can_shade = VARIABLES.CANOPY.Tl_can_shade + 273;

if plants 
%    Tl_net = nanmean(mean([Tl_can_sun  Tl_can_shade]));
%    Tl_net2 = nanmean([nansum(Tl_can_sun.*LAIsun)/sum(LAIsun),nansum(Tl_can_shade.*LAIshade)/sum(LAIshade)]);
    Tl_net_plants = nanmean(nansum([Tl_can_sun.*LAIsun Tl_can_shade.*LAIshade],2)./(LAIshade+LAIsun));
    Tl_net2 = ((Tl_net_plants)*LAI + Tsurf)/(LAI+1);
    
    Tl_net = (LWout/bol/epsv)^(1/4);
else
%    Tl_net = Tsurf;
    Tl_net2 = Tsurf;
    Tl_net = (LWout/bol/epsv)^(1/4);
end
%Ts=Ts+273.15;

%**************************************************************************
%**************************************************************************
%**************************************************************************
%                         I.  FOR EACH LAYER
%**************************************************************************
%**************************************************************************
%**************************************************************************

%**************************************************************************
%                         I.1. Canopy
%**************************************************************************

% INCOME
% Create vector with emissivities
emm=ones(1,2*nlayers)*epsv;
emm(2*nlayers+1)=epsa;
emm(2*nlayers+2)=epss;
emm=emm(:);
% Compute X(e) using the emissivity vector

if (entropymethod == 1)  % Landsberg's Method
    X=0.9652+0.2777*log(1./emm)+0.0511*emm;
elseif (entropymethod == 2)
    X = (1-(c2-c3.*emm).*(45/(4*pi^(4))).*log(emm));
end

% Compute vector of temperatures as column vector [tsun tshade]
vectorT(1:nlayers)=Tl_can_sun;
vectorT(nlayers+1:2*nlayers)=Tl_can_shade;
vectorT(2*nlayers+1)=Tatop;
vectorT(2*nlayers+2)=Tsurf;
vectorT=vectorT(:);

% loop to compute the entropy in each layer
for i=1:1:nlayers
   Ener= LWabs_canM(i,:);
   vec=(4/3).*(Ener').*X./vectorT;   
   ind=~isnan(vec);
   SScanLW_in(i)=sum(vec(ind));     
end


% OUT (EMITTED)

% create vector with emissivities 
emm=ones(nlayers,1)*epsv;

if (entropymethod == 1)  % Landsberg's Method
    X=0.9652+0.2777*log(1./emm)+0.0511*emm;
elseif (entropymethod == 2)
    X = (1-(c2-c3.*emm).*(45/(4*pi^(4))).*log(emm));
end

% Compute entropy and replace nan numbers for zero
SScanLW_out_sun = (4/3).*LWemit_sun./Tl_can_sun.*X*(1-retainLW); 
%SScanLW_out_sun = 2*(4/3).*emm.*bol.*Tl_can_sun.^(3).*LAIsun.*X*(1-retainLW); 
SScanLW_out_sun(isnan(SScanLW_out_sun)) = 0;

SScanLW_out_shade = (4/3).*LWemit_shade./Tl_can_shade.*X*(1-retainLW);
%SScanLW_out_shade = 2*(4/3).*emm.*bol.*Tl_can_shade.^(3).*LAIshade.*X*(1-retainLW); 
SScanLW_out_shade(isnan(SScanLW_out_shade)) = 0;

SScanLW_out = SScanLW_out_sun + SScanLW_out_shade;

%**************************************************************************
%                         I.2. Soil
%**************************************************************************
% Create vector with emissivities
    emm=ones(1,2*nlayers)*epsv;
    emm(2*nlayers+1)=epsa;
    emm(2*nlayers+2)=epss;
    emm=emm(:);
    % Compute X(e) using the emissivity vector

    if (entropymethod == 1)  % Landsberg's Method
        X=0.9652+0.2777*log(1./emm)+0.0511*emm;
    elseif (entropymethod == 2)
        X = (1-(c2-c3.*emm).*(45/(4*pi^(4))).*log(emm));
    end

    % Compute vector of temperatures as column vector [tsun tshade]
    vectorT(1:nlayers)=Tl_can_sun;
    vectorT(nlayers+1:2*nlayers)=Tl_can_shade;
    vectorT(2*nlayers+1)=Tatop;
    vectorT(2*nlayers+2)=Tsurf;
    vectorT=vectorT(:);

    % compute the entropy in each layer
    Ener= LWabs_soilM;
    vec=(4/3).*(Ener').*X./vectorT;   
    ind=~isnan(vec);
    SSsoilLW_in=sum(vec(ind));

    % OUT (EMITTED)
    emm = epss;
    %Esoil = emm*boltz*Tsurf^(4);
    Esoil =LWemit_soil;

    if (entropymethod == 1)  % Landsberg's Method
        X=0.9652+0.2777*log(1./emm)+0.0511*emm;
    elseif (entropymethod == 2)
        X = (1-(c2-c3.*emm).*(45/(4*pi^(4))).*log(emm));
    end

    SSsoilLW_out = (4/3).*Esoil.*X/Tsurf;
%**************************************************************************
%**************************************************************************
%**************************************************************************
%               II.  NET CANOPY
%**************************************************************************
%**************************************************************************
%**************************************************************************
% IN
emm = epsa;

if (entropymethod == 1)  % Landsberg's Method
    X=0.9652+0.2777*log(1./emm)+0.0511*emm;
elseif (entropymethod == 2)
    X = (1-(c2-c3.*emm).*(45/(4*pi^(4))).*log(emm));
end

SSnetLW_Xin = X;
SSnetLW_in = (4/3).*LWin.*X/Tatop;

% OUT (EMITTED)
emm = (epsv*LAI + epsv)/(LAI+1);        % Ponderate between soil and vegetation  

if (entropymethod == 1)  % Landsberg's Method
    X=0.9652+0.2777*log(1./emm)+0.0511*emm;
elseif (entropymethod == 2)
    X = (1-(c2-c3.*emm).*(45/(4*pi^(4))).*log(emm));
end
SSnetLW_Xout = X;
SSnetLW_out = (4/3).*LWout.*X/Tl_net;


%**************************************************************************
%**************************************************************************
%**************************************************************************



% AT EACH LAYER
% canopy
SSresults.SScanLW_in = SScanLW_in;                                         % [W/m^2/K] Entropy due to incoming LW in canopy                                       
SSresults.SScanLW_in_tot = sum(SScanLW_in);                                % [W/m^2/K] Total Entropy due to incoming LW in canopy 
SSresults.SScanLW_out = SScanLW_out;                                       % [W/m^2/K] Entropy due to outgoing LW in canopy
SSresults.SScanLW_out_tot = sum(SScanLW_out);                              % [W/m^2/K] Total Entropy due to incoming LW in canopy 

% soil
SSresults.SSsoilLW_in = SSsoilLW_in;                                       % [W/m^2/K] Entropy due to incoming LW in soils 
SSresults.SSsoilLW_out = SSsoilLW_out;                                     % [W/m^2/K] Entropy due to outgoing LW in soils 

% NET CANOPY
SSresults.SSnetLW_in = SSnetLW_in;                                         % [W/m^2/K] Net (All Ecocystem) Entropy, due to incoming LW   
SSresults.SSnetLW_Xin = SSnetLW_Xin;                                       % [] Net (All Ecocystem) X_LW function for computation of entropy in incoming LW radiation

SSresults.SSnetLW_out = SSnetLW_out;                                       % [W/m^2/K] Net (All Ecocystem) Entropy, due to outgoing LW
SSresults.SSnetLW_Xout = SSnetLW_Xout;                                     % [] Net (All Ecocystem) X_LW function for computation of entropy in outgoing SW radiation                             

SSresults.Tl_net = Tl_net;                                                 % [K] Net (All Ecosystem) Temperature 1 Based on LW
SSresults.Tl_net2 = Tl_net2;                                               % [K] Net (All Ecosystem) Temperature 2 Based on average from canopy and soils 

SSresults.LWin_net = LWin;                                                 % [W/m^2] Total (All Ecosystem) direct incoming SW radiation energy 
SSresults.LWemi_net = sum(LWemit_shade + LWemit_sun)*(1-retainLW) + LWemit_soil; % [W/m^2] Total (All Ecosystem) Emitted LW radiation
SSresults.LWout_net = LWout;                                               % [W/m^2] Total (All Ecosystem) outgoing LW radiation 


