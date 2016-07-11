function [krad,kax,etr] = conductivities(nl_soil,z,dz,zn,rootfr,roottr,...
                                    rpp,smp,vol_liq,eff_porosity,...
                                    K_rad,K_axs,hr,rhc,etr,PLC,kpar_ax,rtcond)
                                
% Compute the radial and axial conductivities of the roots
%-------------------------------------------------------------------------
% Root conductativities (both radial and axial) for each soil layer are
% obtained by weighting the conductivity of the root system by the root
% distribution within the layer. The effect of soil moisture on root
% conductivity is also taken into account.

% Written by Juan Camilo Quijano, UIUC, 2013
%
%------------------------- Input Variables -------------------------------
%       nl_soil     % Number of layers 
%       z           % [mm] Vector with depth of midpoint of each layer 
%       dz          % [mm] Vector with Layer thicknesses information
%       zn          % [mm] Vector with depth of node layer
%       rootfr      % Fine root fraction in each layer
%       roottr      % Tap root fraction in each layer
%       rpp         % [mm] root water potential
%       smp         % [mm] soilwater potential
%       vol_liq     % soil moisture
%       eff_porosity% effective porosity
%       K_rad       % [1/s] Total radial conductivity
%       K_axs       % [mm/s] Total axial conductivity
%       hr          % [1. HR 0. NHR] Logial variable for hr 
%       rhc         % [1. yes 0. no] Root conductivity increases with depth
%       etr         % Transpiration [mm/s]
%       PLC         % [1.Dont work 2. Work] Vectors of 0 or 1 showing layer where radial roots work and does not work      
%       kpar_ax     % Parameter describing the exponential reduction when rtcond = 1
%       rtcond      % [1. Juan's Method, 2. Amenu Method]  Method for the computation of the axial conductivity
%
%------------------------- Output Variables ------------------------------
%       krad        % [/s] radial conductivity of the root 
%       kax         % [mm/s] axial conductivity of the root
%       etr         % Transpiration [mm/s]
%-------------------------------------------------------------------------

    % for radial conductivity, root fraction is used as a weighting factor,
    % because the uptake from a layer is proportional to the total surface 
    % area of roots in that layer. Thus, 

    % for axial conductivity, root density is used as a weighting factor,
    % because the flow is proportional to the total x-sectional area of  
    % the roots in that layer. Thus,                       

% Dereference block
    
    
krad = rootfr.*K_rad.*vol_liq./eff_porosity.*(1-PLC);

if rtcond
    axdist = exp(-kpar_ax*zn);
    %kax = axdist*K_axs./dz;
    kax = K_axs./dz;    
else
    kax = roottr./dz.*K_axs.*vol_liq./eff_porosity;   % Initial Amenu Parameterization
end


% For the case where the root hydraulic conductivity is allowed to increase
% with depth, a linear increasing effect is considered
if rhc == 1     % if conductivity is allowed to increases with depth
    krad = krad + (z/1000).*krad;   % linearly increasing 
    kax  = kax + (z/1000).*kax;
end                               

if hr == 0
   vec=(smp<rpp);
   krad(vec)=0;
   if max(krad)==0
       etr=0;
   end   
end   
   
