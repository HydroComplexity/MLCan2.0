function [vars, km, kl, var] = characteristic(var,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,type)

%=========================================================================
% This function implemenes the characteristic function of Brooks and   
% Corey Model. With the matric pressure the characteristic function is 
% used to computed the soil moisture and hydraulic conductivity.
%
% Written by Juan Camilo Quijano, UIUC, 2008
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

    mc=(var>pentry); %varm(mc)=pentry(mc);   
    nmc=(var<=pentry);
    
    vars(nmc)=thetas(nmc).*(var(nmc)./pentry(nmc)).^(-1./bpar(nmc));
    vars=vars(:);
    kl(nmc)=ks(nmc).*((vars(nmc))./thetas(nmc)).^(2.*bpar(nmc)+3);  
    
    vars(mc) = thetas(mc);
    kl(mc) = ks(mc);
    vars=vars(:);
    
    [varsP] = comave(vars,zsoi,dzsoi,zisoi,nl_soil);[bparP] = comave(bpar,zsoi,dzsoi,zisoi,nl_soil);  
    [thetasP] = comave(thetas,zsoi,dzsoi,zisoi,nl_soil);[ksP] = comave(ks,zsoi,dzsoi,zisoi,nl_soil);
    km=ksP.*((varsP)./thetasP).^(2.*bparP+3);   

elseif type ==2 % compute psi based on theta
%   check if psi is smaller than saturated    
%   mc=(var>thetas); var(mc)=thetas(mc);          
    vars=pentry.*(var./thetas).^(-bpar);
    kl=ks.*((var)./thetas).^(2.*bpar+3);                

    [varP] = comave(var,zsoi,dzsoi,zisoi,nl_soil);[bparP] = comave(bpar,zsoi,dzsoi,zisoi,nl_soil);  
    [thetasP] = comave(thetas,zsoi,dzsoi,zisoi,nl_soil);[ksP] = comave(ks,zsoi,dzsoi,zisoi,nl_soil);
    km=ksP.*((varP)./thetasP).^(2.*bparP+3);            
end
        
        
        