function P=PEDO_TRANSFER_FUNCTIONS(S,C,OM)

% Pedotransfer function from Saxon and Rawls 2006
% Soil Water Characteristic Estimates by Texture and Organic Matter for
% Hydrologic Solutions, Soil Sci. Soc. Am. J. 70:1569–1578 (2006).
% doi:10.2136/sssaj2005.0117

% Inputs:
% S  % Sand by mass
% C  % Clay by mass
% OM % Organic matter by mass


Theta_1500t = -0.024.*S + 0.487.*C + 0.006.*OM ...
    + 0.005.*(S.*OM) - 0.013.*(C.*OM) ...
    + 0.068.*(S.*C) + 0.031;

Theta_1500 = Theta_1500t + (0.14.*Theta_1500t - 0.02);

Theta_33t = -0.251.*S + 0.195.*C + 0.011.*OM ...
    + 0.006.*(S.*OM) - 0.027.*(C.*OM) ...
    + 0.452.*(S.*C) + 0.299;

Theta_33 = Theta_33t + (1.238.*Theta_33t.*Theta_33t - 0.374.*Theta_33t - 0.015);

Theta_Sm33t=0.278.*S + 0.034.*C + 0.022.*OM...
    -0.018.*(S.*OM)-0.027.*(C.*OM)...
    -0.584.*(S.*C) + 0.078;...
    
Theta_Sm33=Theta_Sm33t+(0.636.*Theta_Sm33t-0.107);

Theta_sat=Theta_33+Theta_Sm33-0.097.*S+0.043;


Phi_et=-21.67*S-27.93*C-81.97*Theta_Sm33...
    +71.12*(S.*Theta_Sm33)+8.29*(C.*Theta_Sm33)...
    +14.05*(S.*C)+27.16;

Phi_e=Phi_et+(0.02*Phi_et.^2-0.113*Phi_et-0.7);

B=(log(1500)-log(33))./(log(Theta_33)-log(Theta_1500));

lam=1./B;
Ksat=1930.*(Theta_sat-Theta_33).^(3-lam); %mm/hr

BCbeta=1./(3+2./lam);



P.Ksat	 =Ksat;				% Saturated conductivity (mm/hr)
P.Sr	 =0;				% Residual VWC, with respect to percolation behaviour: K(VWCr)=0
P.Ss	 =Theta_sat;		% Saturation water content i.e. porosity
P.BCbeta	=BCbeta;        % Brooks Corey Shape Parameter
P.lam = lam;                % Clapp-Hornberger "b" parameter [-]
P.hb = - Phi_e/1000;        % Air entry pressure in MPa
P.S_init=Theta_sat/2;		% Initial VWC
P.S_fc=Theta_33;			% Storage above which evaporation is at potential rate
P.S_wp=Theta_1500;			% Storage below which evaporation is zero
