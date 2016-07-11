function [vars, km, kl, var] = characteristic_vanGenuchten(var,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,type,alpha,n,lambda,m,thetar)

%=========================================================================
% This function implemenes the characteristic function of VanGenuchten
% With the matric pressure the characteristic function is 
% used to computed the soil moisture and hydraulic conductivity.
%
% Written by Dongkook Woo, UIUC, 2016
%
%------------------------- Input Variables -------------------------------
%       var             % Incoming variable. (Matric Pressure/Soil Moisture)
%       thetas          % saturated soil moisture 
%       ks              % Hydraulic conductivity at saturation [m/s]
%       nl_soil         % Number of layers in the soil        
%       zsoi            % Vector with depth of midpoint of each layer 
%       dzsoi           % Vector with Layer thicknesses information
%       zisoi           % Vector with depth of bounries of layers
%       bpar            % Clapp-Hornbereger "b" parameter
%       pentry          % [mm] minimum soil suction, i.e., soil potential at saturation 
%       type            % [1] Compute theta based on psi, [2] compute psi
%                          based on theta
%
%------------------------- Output Variables ------------------------------
%       vars            % Outgoing variable. (Matric Pressure/Soil Moisture)
%       km              % Hydraulic conductivity for interface of each layer
%       kl              % Hydraulic conductivity for each layer 
%       var             % Incoming variable  
%-------------------------- Unit Conversions -----------------------------

%CASE BROOKS AND COREY TYPE II
% Constants for Brooks and Corey
if type==1; % compute theta based on psi
%   check if theta is smaller than saturated

    vars = (thetas - thetar)./(1 + (alpha.*abs(var)).^n).^m + thetar;
    vars(var>=0.0) = thetas(var>=0.0);
    
    Se = ((vars - thetar)./(thetas - thetar));
    kl = ks.*Se.^(1/2).*(1 - (1 - Se.^(1./m)).^m).^2;
    
    [SeP] = comave(Se,zsoi,dzsoi,zisoi,nl_soil); [mP] = comave(m,zsoi,dzsoi,zisoi,nl_soil);  
    [ksP] = comave(ks,zsoi,dzsoi,zisoi,nl_soil);
    km=ksP.*SeP.^(1/2).*(1 - (1 - SeP.^(1./mP)).^mP).^2;
    
elseif type ==2 % compute psi based on theta
%    theta = (thetas - thetar)./(1 + (alpha.*abs(h)).^n).^m + thetar;
%    var = (thetas - thetar)./(1 + (alpha.*abs(vars)).^n).^m + thetar;
    
    vars = -(((1./((var-thetar)./(thetas-thetar)).^(1./m))-1).^(1./n))./alpha;
    
    Se = ((var - thetar)./(thetas - thetar));
    
    % Compute the hydraulic conductivity [eqn 8]
    kl = ks.*Se.^(1/2).*(1 - (1 - Se.^(1./m)).^m).^2;
    
    [SeP] = comave(Se,zsoi,dzsoi,zisoi,nl_soil); [mP] = comave(m,zsoi,dzsoi,zisoi,nl_soil);  
    [ksP] = comave(ks,zsoi,dzsoi,zisoi,nl_soil);
    km=ksP.*SeP.^(1/2).*(1 - (1 - SeP.^(1./mP)).^mP).^2;
end