
function [SSresults] = ENTROPY_SW (SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
                                          SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
                                          SWout, fdiff, entropymethod, Rg, zicesl, PARAMS)

%=========================================================================
% This code computes the fluxes of entropy due to shortwave. 
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       SWcandir_in                    % [W/m2] Direct Incoming Shortwave Radiation in canopy
%       SWcandir_out                   % [W/m2] Direct Outgoing Shortwave Radiation in canopy
%       SWcandif_in                    % [W/m2] Diffuse Incoming Shortwave Radiation in canopy 
%       SWcandif_out                   % [W/m2] Diffuse Outgoing Shortwave Radiation in canopy
%       SWsoildir_in                   % [W/m2] Direct Incoming Shortwave Radiation in soil
%       SWsoildir_out                  % [W/m2] Direct Outgoing Shortwave Radiation in soil
%       SWsoildif_in                   % [W/m2] Diffuse Incoming Shortwave Radiation in soil
%       SWsoildif_out                  % [W/m2] Diffuse Outgoing Shortwave Radiation in soil
%       SWout                          % [W/m2] Total Outgoing Shortwave Radiation 
%       fdiff                          % [] Canopy Fraction that receives diffuse
%       entropymethod                  % [] Method used to compute the entropy in radiation
%       Rg                             % [W/m2] Incoming SW radiation 
%       zicesl                         % [mm] Ice water depth in Snow-litter pack 
%       PARAMS                         % PARAMS structure
%------------------------- Output Variables ------------------------------
%       SSresults                      % SSresults structure 
% 
%========================================================================
% De reference block
%*********************************************************************
                                                                            
% PARAMETERS TO COMPUTE GREYBODY ENTROPY (Wright et al 2001. Int J of Eng.)
c2 = PARAMS.Entropy.c2;                 %[] Constant to conpute entropy SW 
c3 = PARAMS.Entropy.c3;                 %[] Constant to conpute entropy SW  
nlayers=PARAMS.CanStruc.nl_can;         %[] Number of layers in soil structure
                                      
To=5760 ;                               %[K] Solar Temperatute 

%**************************************************************************
%**************************************************************************
%**************************************************************************
%                        I.  FOR EACH LAYER
%**************************************************************************
%**************************************************************************
%**************************************************************************

%**************************************************************************
%                       I.1. canopy
%**************************************************************************
% Computation of entropy due to direct shortwave coming
k= 2.31*10^(-4);    %[s1/e1]  [1/�K]
                    %        s1:Solar Constant of second kind 
                    %        e1:Solar Constant of first kind
SScandir_in = (SWcandir_in) * k;
SScandir_out = (SWcandir_out) * k;

%**************************************************************************
% Computation of entropy due to diffuse shortwave coming 
% Diffuse in 
kin=SWcandif_in./pi;
kout=SWcandif_out./pi; 
ko= 1.99*10^(7); % Solar Energy Radiation in extraterrestrial space [J/s/m2]
eein=kin./ko;
eeout=kout./ko;
    
if (Rg>0 && ~isempty(fdiff) && entropymethod == 1)  % Landsberg's Method

    Xin=0.9652+0.2777*log(1./eein)+0.0511*eein;            
    SScandif_in=(4/3).*(SWcandif_in/To).*Xin; 
    SScandif_in(isnan(SScandif_in)) = 0;        % Correct for nan numbers from eein=0
    
    Xout=0.9652+0.2777*log(1./eeout)+0.0511*eeout;            
    SScandif_out=(4/3).*(SWcandif_in/To).*Xout;        
    SScandif_out(isnan(SScandif_out)) = 0;        % Correct for nan numbers from eein=0

elseif (Rg>0 && ~isempty(fdiff) && entropymethod == 2) % Wright's Method

    Xin =  (1-(c2-c3*eein).*(45./(4*pi^(4))).*log(eein));
    SScandif_in=(4/3).*(SWcandif_in/To).*Xin;        
    SScandif_in(isnan(SScandif_in)) = 0;        % Correct for nan numbers from eein=0

    Xout = (1-(c2-c3*eeout).*(45./(4*pi^(4))).*log(eeout));
    SScandif_out=(4/3).*(SWcandif_in/To).*Xout;        
    SScandif_out(isnan(SScandif_out)) = 0;        % Correct for nan numbers from eein=0

else
    
    SScandif_in = zeros(nlayers,1);
    SScandif_out = zeros(nlayers,1);

end

%**************************************************************************
%                       I.2. soil
%**************************************************************************
% Computation of entropy due to direct shortwave coming
    k= 2.31*10^(-4);    %[s1/e1]  [1/�K]
                        %        s1:Solar Constant of second kind 
                        %        e1:Solar Constant of first kind
    SSsoildir_in = (SWsoildir_in) * k;
    SSsoildir_out = (SWsoildir_out) * k;
    
%**************************************************************************
% Computation of entropy due to diffuse shortwave coming 
% Diffuse in 
    kin=SWsoildif_in./pi;
    kout=SWsoildif_out./pi; 
    ko= 1.99*10^(7); % Solar Energy Radiation in extraterrestrial space [J/s/m2]
    eein=kin./ko;
    eeout=kout./ko;

    if (Rg>0 && ~isempty(fdiff) && entropymethod == 1)  % Landsberg's Method

        Xin=0.9652+0.2777*log(1./eein)+0.0511*eein;            
        SSsoildif_in=(4/3).*(SWsoildif_in/To).*Xin;  
        SSsoildif_in(isnan(SSsoildif_in)) = 0;        % Correct for nan numbers from eein=0    

        Xout=0.9652+0.2777*log(1./eeout)+0.0511*eeout;            
        SSsoildif_out=(4/3).*(SWsoildif_in/To).*Xout;
        SSsoildif_out(isnan(SSsoildif_out)) = 0;        % Correct for nan numbers from eein=0


    elseif (Rg>0 && ~isempty(fdiff) && entropymethod == 2) % Wright's Method

        Xin = (1-(c2-c3*eein)*(45/(4*pi^(4)))*log(eein));
        SSsoildif_in=(4/3).*(SWsoildif_in/To).*Xin;        
        SSsoildif_in(isnan(SSsoildif_in)) = 0;        % Correct for nan numbers from eein=0    

        Xout = (1-(c2-c3*eeout)*(45/(4*pi^(4)))*log(eeout));
        SSsoildif_out=(4/3).*(SWsoildif_in/To).*Xout;        
        SSsoildif_out(isnan(SSsoildif_out)) = 0;        % Correct for nan numbers from eein=0

    else

        SSsoildif_in = 0;
        SSsoildif_out = 0;

    end
%**************************************************************************
%**************************************************************************
%**************************************************************************
%                         II.  THE NET CANOPY
%**************************************************************************
%**************************************************************************
%**************************************************************************

SWnetdir_in = Rg*(1-fdiff);
SWnetdif_in = Rg*(fdiff);

SWnetdir_out = 0;
SWnetdif_out = SWout;

%**************************************************************************
% Computation of entropy due to direct shortwave coming
k= 2.31*10^(-4);    %[s1/e1]  [1/�K]
                    %        s1:Solar Constant of second kind 
                    %        e1:Solar Constant of first kind
SSnetdir_in = (SWnetdir_in) * k;
SSnetdir_out = (SWnetdir_out) * k;

%**************************************************************************
% Computation of entropy due to diffuse shortwave coming 
% Diffuse in 

kin=SWnetdif_in./pi;
kout=SWnetdif_out./pi;  % Solar Energy Radiation in extraterrestrial space [J/s/m2]
ko= 1.99*10^(7);        % Solar Energy Radiation in extraterrestrial space [J/s/m2]
eein=kin./ko;
eeout=kout./ko;
    
if (Rg>0 && ~isempty(fdiff) && entropymethod == 1)  % Landsberg's Method

    Xin=0.9652+0.2777*log(1./eein)+0.0511*eein;            
    SSnetdif_in=(4/3).*(SWnetdif_in/To).*Xin;
    SSnetdif_in(eein == 0) = 0;
    
    Xout=0.9652+0.2777*log(1./eeout)+0.0511*eeout;            
    SSnetdif_out=(4/3).*(SWnetdif_in/To).*Xout;        
    SSnetdif_out(eeout == 0) = 0;
    
elseif (Rg>0 && ~isempty(fdiff) && entropymethod == 2) % Wright's Method

    Xin = (1-(c2-c3*eein)*(45/(4*pi^(4)))*log(eein));
    SSnetdif_in=(4/3).*(SWnetdif_in/To).*Xin;
    SSnetdif_in(eein == 0) = 0;    
    
    Xout = (1-(c2-c3*eeout)*(45/(4*pi^(4)))*log(eeout));
    SSnetdif_out=(4/3).*(SWnetdif_in/To).*Xout;        
    SSnetdif_out(eeout == 0) = 0;

else    
    SSnetdif_in = 0;
    SSnetdif_out = 0;
    Xin = nan;
    Xout = nan;
end
SSnetSW_Xin = Xin;
SSnetSW_Xout = Xout;


%**************************************************************************
%**************************************************************************
%**************************************************************************
%               III.  STORE ENTROPY CALCULATIONS IN STRUCTURE
%**************************************************************************
%**************************************************************************
%**************************************************************************

SSresults.SScandir_in = SScandir_in;                                        % [W/m^2/K] Entropy due to incoming SW direct in canopy
SSresults.SScandir_in_tot = sum(SScandir_in);                               % [W/m^2/K] Total Entropy due to incoming SW direct in canopy
SSresults.SScandir_out = SScandir_out;                                      % [W/m^2/K] Entropy due to outgoing SW direct in canopy    
SSresults.SScandir_out_tot = sum(SScandir_out);                             % [W/m^2/K] Total Entropy due to outgoing SW direct in canopy 
SSresults.SScandif_in = SScandif_in;                                        % [W/m^2/K] Entropy due to incoming SW diffuse in canopy
SSresults.SScandif_in_tot = sum(SScandif_in);                               % [W/m^2/K] Total Entropy due to incoming SW diffuse in canopy
SSresults.SScandif_out = SScandif_out;                                      % [W/m^2/K] Entropy due to outgoing  SW diffuse in canopy
SSresults.SScandif_out_tot = sum(SScandif_out);                             % [W/m^2/K] Total Entropy due to outgoing  SW diffuse in canopy

SSresults.SSsoildir_in = SSsoildir_in;                                      % [W/m^2/K] Entropy due to incoming  SW direct in soil
SSresults.SSsoildir_out = SSsoildir_out;                                    % [W/m^2/K] Entropy due to outgoing  SW direct in soil    
SSresults.SSsoildif_in = SSsoildif_in;                                      % [W/m^2/K] Entropy due to incoming  SW diffuse in soil
SSresults.SSsoildif_out = SSsoildif_out;                                    % [W/m^2/K] Entropy due to outgoing  SW diffuse in soil

SSresults.SSnetdir_in = SSnetdir_in;                                        % [W/m^2/K] Net (All Ecocystem) Entropy, due to incoming SW direct
SSresults.SSnetdir_out = SSnetdir_out;                                      % [W/m^2/K] Net (All Ecocystem) Entropy, due to outgoing SW direct
SSresults.SSnetdif_in = SSnetdif_in;                                        % [W/m^2/K] Net (All Ecocystem) Entropy, due to incoming SW diffuse
SSresults.SSnetdif_out = SSnetdif_out;                                      % [W/m^2/K] Net (All Ecocystem) Entropy, due to outgoing SW diffuse

SSresults.SWdir_in = SWnetdir_in;                                           % [W/m^2] Total direct incoming SW radiation energy 
SSresults.SWdir_out = 0;                                                    % [W/m^2] Total direct outgoing SW radiation energy 
SSresults.SWdif_in = SWnetdif_in;                                           % [W/m^2] Total diffuse incoming SW radiation energy 
SSresults.SWdif_out = SWout;                                                % [W/m^2] Total diffuse outgoing SW radiation energy 

SSresults.SSnetSW_Xin = SSnetSW_Xin;                                        % [] Net (All Ecocystem) X_SW function for computation of entropy in incoming SW radiation
SSresults.SSnetSW_Xout = SSnetSW_Xout;                                      % [] Net (All Ecocystem) X_SW function for computation of entropy in outgoing SW radiation



