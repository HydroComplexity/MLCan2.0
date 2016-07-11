function [C] = computecS_vanGenuchten(vpsi,vtheta,thetas,nl_soil,bpar,pentry,alpha,n,lambda,m,thetar)                       
%=========================================================================
% This code computes C which is defined as the rate of change of soil
% moisture with respect to soil matric pressure. The analytical solution
% of the derivative of Brooks and Corey equation is implemented.
%
% Written by Dongkook Woo, UIUC, 2016
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

% Compute the specific moisture storage (derivative of eqn 21:
% C = d(theta)/dh
C = -alpha.*n.*sign(vpsi).*(1./n - 1).*(alpha.*abs(vpsi)).^(n - 1).*(thetar - thetas).*((alpha.*abs(vpsi)).^n + 1).^(1./n - 2);
C(vpsi>=0.0) = 0;