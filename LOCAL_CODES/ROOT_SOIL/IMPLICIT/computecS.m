function [C] = computecS(vpsi,vtheta,thetas,nl_soil,bpar,pentry)                       
%=========================================================================
% This code computes C which is defined as the rate of change of soil
% moisture with respect to soil matric pressure. The analytical solution
% of the derivative of Brooks and Corey equation is implemented.
%
% Written by Juan Camilo Quijano, UIUC, 2008
%
%------------------------- Input Variables -------------------------------
%       vpsi     % [mm] soil matric pressure
%       vpsi     % [mm] soil water content
%       thetas   % [] saturated soil moisture
%       nl_soil  % [] number of layers  
%       bpar     % Clapp-Hornbereger "b" parameter
%       pentry   % [mm] minimum soil suction, i.e., soil potential at
%                   saturation 
%------------------------- Output Variables ------------------------------
%------------------------- Output Variables ------------------------------
%       C        % variable C. Rate of change of soil moisture with respect
%                  to soil matric pressure.
%-------------------------------------------------------------------------   

%mc=vpsi>pentry;vpsi(mc)=pentry(mc);
psicr = pentry - 0.01*abs(pentry);

% Lower than critical (power law) 
ind1 = vpsi <= pentry;
C(ind1)=thetas(ind1)./pentry(ind1).*(-1./bpar(ind1)).*(vpsi(ind1)./pentry(ind1)).^(-1./bpar(ind1)-1);

% Higher than critical lower than zero (parabola)
ind2a = vpsi > pentry;
ind2b = vpsi <= 0;
ind2 = logical(ind2a.*ind2b);
thetacr=thetas.*(psicr./pentry).^(-1./bpar);
Wcr = thetacr./thetas; 
m = (psicr)./((1-Wcr).^(2)) - (psicr.*bpar)./(Wcr.*(1-Wcr)); 
n = 2*Wcr - ((psicr.*bpar)./(m.*Wcr)) - 1;
C(ind2) = 1 ./ (-m(ind2).*(((2*vtheta(ind2))./(thetas(ind2).^2))   -    1./(thetas(ind2))   -  n(ind2)./thetas(ind2) ));


% Higher than zero
ind3 = vpsi > 0;
C(ind3) = 0;
%