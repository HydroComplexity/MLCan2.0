
function [SSresults] = ENTROPY_PHO (VARIABLES, SWITCHES, VERTSTRUC, PARAMS, SSresults)



%=========================================================================
% This code computes the fluxes of entropy due to photosynthesis. Onlyt 
% the incoming fluxes are computed. It is assumed that all the energy
% captured through photosynthesis is released as heat
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES                       % VARIABLES structure
%       SWITCHES                        % SWITCHES structure
%       VERTSTRUC                       % SWITCHES structure
%       PARAMS                          % SWITCHES structure
%       SSresults                       % SSresults structure 
%------------------------- Output Variables ------------------------------
%       SSresults                      % SSresults structure 
% 
%========================================================================
% De reference block
%*********************************************************************

% SWITCHES
entropymethod = SWITCHES.entropymethod;                                    % [] Method used to compute the entropy in radiation

% PARAMS
c2 = PARAMS.Entropy.c2;                                                    %[] Constant to conpute entropy SW                         
c3 = PARAMS.Entropy.c3;                                                    %[] Constant to conpute entropy SW 


To=5760 ; % Solar Temperatute [K]

%**************************************************************************
% Photosynthetic fluxes
%**************************************************************************
LAI_sun = VARIABLES.CANOPY.LAIsun;
LAI_shade = VARIABLES.CANOPY.LAIshade;
fLAIz = VERTSTRUC.fLAIz;
An_shade = VARIABLES.CANOPY.An_shade;
An_sun = VARIABLES.CANOPY.An_sun;

An_shade_prof = (sum(An_shade.*fLAIz,2)).*LAI_shade;             
An_sun_prof = (sum(An_sun.*fLAIz,2)).*LAI_sun;                   

Epho_dif = max((sum(An_shade_prof))*.506,0);
Epho_dir = max((sum(An_sun_prof))*.506,0);
Epho_tot = (sum(An_sun_prof)+sum(An_shade_prof))*0.506;


% Computation of entropy due to direct shortwave coming into the CV
k= 2.31*10^(-4);    %[s1/e1]  [K]
                    %        s1:Solar Constant of second kind 
                    %        e1:Solar Constant of first kind
SSphodir_in = (Epho_dir) * k;

% Computation of entropy due to diffuse shortwave coming into the CV

% Diffuse in 
kin=Epho_dif./pi;
ko= 1.99*10^(7); % Solar Energy Radiation in extraterrestrial space [J/s/m2]
eein=kin./ko;
    
if entropymethod == 1 % Landsberg's Method

    Xin=0.9652+0.2777*log(1./eein)+0.0511*eein;            
    SSphodif_in=(4/3).*(Epho_dif/To).*Xin; 
    SSphodif_in(isnan(SSphodif_in)) = 0;        % Correct for nan numbers from eein=0
    
elseif entropymethod == 2 % Wright's Method

    Xin =  (1-(c2-c3*eein).*(45./(4*pi^(4))).*log(eein));
    SSphodif_in=(4/3).*(Epho_dif/To).*Xin;        
    SSphodif_in(isnan(SSphodif_in)) = 0;        % Correct for nan numbers from eein=0

end
if Epho_dif>0
    stop = 213;
end

% Store the variables 
SSresults.Epho_dif = Epho_dif;                                             % [W/m2] Energy of difussed SW being captured through Photosynthesis 
SSresults.Epho_dir = Epho_dir;                                             % [W/m2] Energy of direct SW being captured through Photosynthesis
SSresults.Epho_tot = Epho_tot;                                             % [W/m2] Total Energy of direct SW being captured through Photosynthesis
SSresults.SSphodir_in = SSphodir_in;                                       % [W/m2/K] Entropy flux fue to direct SW capture through Photosynthesis
SSresults.SSphodif_in = SSphodif_in;                                       % [W/m2/K] Entropy flux fue to direct SW capture through Photosynthesis


