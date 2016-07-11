function [VERTSTRUC] = SOIL_PROPERTIES(PARAMS, VERTSTRUC)

%=========================================================================
% This function calculates the profiles of the soil properties that are
% dependent solely on soil composition (sand and clay content)
%
% See Oleson et al (2004) (CLM Documentation) for equation # references
%
% Written by Darren Drewry
%=========================================================================

% Dereference Structure Values
    nl_soil = VERTSTRUC.nl_soil;

    sand = VERTSTRUC.sand*ones(nl_soil,1);
    clay = VERTSTRUC.clay*ones(nl_soil,1);
%    load sand.mat;
%    load clay.mat;
    
    zhs = VERTSTRUC.zhs;
    smpmin = PARAMS.Soil.smpmin;
    scalek = PARAMS.Soil.scalek;

% HEAT CAPACITY OF SOIL SOLIDS [J / m^3 / K]
    HC_sol = (2.128*sand+2.385*clay) ./ (sand+clay) * 10^6;     % (6.68) 
    
% POROSITY = WATER CONTENT AT SATURATION [-]
    porsl = 0.489 - 0.00126*sand;                               % (7.72)
    porsl = 0.4 * ones(size(sand));
     
% MINIMUM SOIL SUCTION = SOIL POTENTIAL AT SATURATION [mm]
    psi0 = -10 * ( 10.^(1.88-0.0131*sand) );                    % (7.75)
     
% Clapp-Hornberger "b" parameter [-]
    bsw = 2.91 + 0.159*clay;                                    % (7.73)
    
% THERMAL CONDUCTIVITY OF SOIL MINERALS [W / m / K]
         TK_sol  = ((8.80*sand+2.92*clay) ./ (sand+clay));             % (6.61)
     
% BULK DENSITY [kg / m^3]
    rhod = 2700*(1 - porsl);                                    % (before 6.62)
    
% THERMAL CONDUCTIVITY OF DRY SOIL [W / m / K]
    TK_dry = (0.135*rhod + 64.7) ./ (2700 - 0.947*rhod);        % (6.62)
     
% HYDRAULIC CONDUCTIVITY AT SATURATION [mm / s] 
    HKsat = 0.0070556 * ( 10.^(-0.884+0.0153*sand) )...         % (7.71)
                      .*exp(-zhs*scalek);
%     HKsat = 0.0070556 * ( 10.^(-0.884+0.0153*sand));         % (7.71)

% SOIL MOISTURE CONTENT AT COMPLETE DRYNESS
    theta_dry = porsl .* ( psi0./smpmin ).^(1./bsw);            % 
  
% Dongkook: Start
% Soil moisture parameter for vanGenuchten
% Minasny and McBratney, 2007

% For n
%S1=1./(1+exp(-(24.547 - 0.238.*sand - 0.082.*clay)));
%S2=1./(1+exp(-(-3.569 + 0.081.*sand)));
%S3=1./(1+exp(-(0.694 - 0.024.*sand + 0.048.*clay)));

%n=2.18+0.11.*(48.087 - 44.954.*S1 - 1.023.*S2-3.896.*S3);

% For alpha
% parameter related to the inverse of the air entry suction [1/cm] to [1/mm]
% Markus Tuller and Dani Or, 2003
% http://www.mathworks.com/matlabcentral/fileexchange/45468-soil-classification/content/soil_classification.m
sandp=sand./100;
clayp=clay./100;

siltp=1-sandp-clayp;
for i=1:length(clayp)
    if (siltp(i)+1.5*clayp(i))<.15
        %SC{i,:}='SAND';
        alpha(i) = 0.035*(1/10);
        n(i) = 3.19;
        thetar(i) = 0.058;
    elseif ((siltp(i)+1.5*clayp(i))>=.15)&((siltp(i)+2*clayp(i))<.3)
        %SC{i,:}='LOAMY SAND';
        alpha(i) = 0.035*(1/10);
        n(i) = 2.39;
        thetar(i) = 0.074;
    elseif (clayp(i)>=0.07) & (clayp(i)<=0.2) & (sandp(i)>0.52) & ((siltp(i)+2*clayp(i))>=0.3)
        %SC{i,:}='SANDY LOAM';
        alpha(i) = 0.021*(1/10);
        n(i) = 1.61;
        thetar(i) = 0.067;
    elseif (clayp(i)<0.07) & (siltp(i) < 0.5) & ((siltp(i)+2*clayp(i))>=0.3)
        %SC{i,:}='SANDY LOAM';
        alpha(i) = 0.021*(1/10);
        n(i) = 1.61;
        thetar(i) = 0.067;
    elseif (clayp(i)>=0.07) & (clayp(i)<=0.27) & (siltp(i)>=0.28) & (siltp(i)<0.5) & (sandp(i)<=0.52)
        %SC{i,:}='LOAM';
        alpha(i) = 0.036*(1/10);
        n(i) = 1.31;
        thetar(i) = 0.083;
    elseif ((siltp(i)>0.5) & clayp(i)>0.12 & clayp(i)<0.27) | (siltp(i)>=0.5 & siltp(i)<0.8 & clayp(i)<0.12)
        %SC{i,:}='SILT LOAM';
        alpha(i) = 0.025*(1/10);
        n(i) = 1.39;
        thetar(i) = 0.061;
    elseif siltp(i)>=0.8 & clayp(i)<0.12
        %SC{i,:}='SILT';
        alpha(i) = 0.006*(1/10);
        n(i) = 1.53;
        thetar(i) = 0.123;
    elseif clayp(i)>=0.2 & clayp(i)<0.35 & siltp(i)<0.28 & sandp(i)>0.45
        %SC{i,:}='SANDY CLAY LOAM';
        alpha(i) = 0.033*(1/10);
        n(i) = 1.49;
        thetar(i) = 0.086;
    elseif clayp(i)>=0.27 & clayp(i) <0.4 & sandp(i)>0.2 & sandp(i)<=0.45
        %SC{i,:}='CLAY LOAM';
        alpha(i) = 0.030*(1/10);
        n(i) = 1.37;
        thetar(i) = 0.129;
    elseif clayp(i)>=0.27 & clayp(i)<0.4 & sandp(i)<=0.2
        %SC{i,:}='SILTY CLAY LOAM';
        alpha(i) = 0.027*(1/10);
        n(i) = 1.41;
        thetar(i) = 0.098;
    elseif clayp(i)>=0.35 & sandp(i)>45
        %SC{i,:}='SANDY CLAY';
        alpha(i) = 0.025*(1/10);
        n(i) = 1.4;
        thetar(i) = 0.08;
    elseif clayp(i)>=0.4 & siltp(i)>=0.4
        %SC{i,:}='SILTY CLAY';
        alpha(i) = 0.023*(1/10);
        n(i) = 1.39;
        thetar(i) = 0.163;
    elseif clayp(i)>= 0.4 & sandp(i)<=0.45 & siltp(i)<0.4
        %SC{i,:}='CLAY'
        alpha(i) = 0.021*(1/10);
        n(i) = 1.20;
        thetar(i) = 0.102;
    end
end
alpha=alpha(:);
n=n(:);
thetar=thetar(:);

%alpha(:)   = 0.01*(1/10);         % parameter related to the inverse of the air entry suction [1/cm] to [1/mm]
%n(:)       = 2.0;                  % Pore-size distributions [-]
% Dongkook: End

%*************************************************************************
% STORE IN STRUCTURE
%*************************************************************************

    % ASSIGN
        VERTSTRUC.HC_sol = HC_sol;
        VERTSTRUC.porsl = porsl;
        VERTSTRUC.psi0 = psi0;
        VERTSTRUC.bsw = bsw;
        VERTSTRUC.TK_sol = TK_sol;
        VERTSTRUC.TK_dry = TK_dry;
        VERTSTRUC.HKsat = HKsat;
        VERTSTRUC.theta_dry = theta_dry;        
        VERTSTRUC.eff_poros = VERTSTRUC.porsl;
        % Dongkook Woo - Edit 
        VERTSTRUC.rhod= rhod;
        VERTSTRUC.clayColum = clay;
        
        VERTSTRUC.VanGen_n=n;
        VERTSTRUC.VanGen_alpha=alpha;
        VERTSTRUC.VanGen_thetar=thetar;
        % Dongkook Woo - Edit End
    
    
    