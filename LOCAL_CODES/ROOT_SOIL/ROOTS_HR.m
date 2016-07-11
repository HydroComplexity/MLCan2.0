function [rpp, rpp_weight, krad, kax] = ROOTS_HR( SWITCHES, VERTSTRUC, PARAMS, VARIABLES )

%=========================================================================
% This code solves the model for water flow in the plant root system. 
% The upper boundary condition is set to the transpiration rate while 
% the lower boundary is set to no flux.
%
%------------------------- Input Variables -------------------------------
%       dtime           % time step [s]
%       z               % layer depth [mm]
%       dz              % layer thickness [mm]
%       zi              % interface level below a "z" level [mm]
%       TR             % actual transpiration [mm/s]
%       rootfr          % root fraction of a layer [-]
%       rpp             % root pressure potential [mm]
%       smp             % soil matric potential [mm]
%       vol_liq         % soil water per unit volume [mm/mm]
%       eff_porosity    % effective porosity of soil
%       thetadry        % soil moisture content at dryness [-]
%       K_rad           % radial conductivity of the root system [s^-1]
%       K_axs           % axial specific conductivity of the root system [mm/s]
%       hr              % option for hydraulic redistribution [-]
%       rhc             % option for root hydraulic conductivity [-]
%
%------------------------- Output Variables ------------------------------
%       rpp             % root pressure potential [mm]
%       krad            % radial conductivity of the root [/s]
%       kax             % axial conductivity of the root [mm/s]
%                    
%-------------------------- local variables ------------------------------
%       amx             % "a" left off diagonal of tridiagonal matrix
%       bmx             % "b" diagonal column for tridiagonal matrix
%       cmx             % "c" right off diagonal tridiagonal matrix
%       dmx             % "d" forcing term of tridiagonal matrix
% 
%=========================================================================

    
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    rhc = SWITCHES.rhc;

    TR = VARIABLES.CANOPY.TR_can;
    
    smp = VARIABLES.SOIL.smp;
    vol_liq = VARIABLES.SOIL.volliq;
    
    z = VERTSTRUC.znsmm;
    dz = VERTSTRUC.dzsmm;
    rootfr = VERTSTRUC.rootfr;
    thetadry = VERTSTRUC.theta_dry;
    eff_porosity = VERTSTRUC.eff_poros;

    K_rad = PARAMS.Soil.K_rad;
    K_axs = PARAMS.Soil.K_axs;
%*************************************************************************
%*************************************************************************
    
rinds = find(rootfr>0);
nrinds = find(rootfr<=0);

N = length(rinds);  % # soil layers

% Compute the radial and axial conductivities of the roots
%-------------------------------------------------------------------------
% Root conductativities (both radial and axial) for each soil layer are
% obtained by weighting the conductivity of the root system by the root
% distribution within the layer. The effect of soil moisture on root
% conductivity is also taken into account.

binds = [];
ginds = [];
for i = 1:N

    if (vol_liq(i) <= thetadry(i) )
        krad(i) = 0;
        binds = [binds, i];
    else
        % for radial conductivity, root fraction is used as a weighting factor,
        % because the uptake from a layer is proportional to the total surface 
        % area of roots in that layer. Thus, 
        krad(i) = rootfr(i)*K_rad ...
                  * vol_liq(i)/eff_porosity(i);
    end
        % for axial conductivity, root density is used as a weighting factor,
        % because the flow is proportional to the total x-sectional area of  
        % the roots in that layer. Thus,                       
        kax(i) = (rootfr(i)/dz(i))*K_axs ...
                 * vol_liq(i)/eff_porosity(i); 
             
        ginds = [ginds, i];
    

end

krad(nrinds) = 0;
kax(nrinds) = 0;

krad = krad(:);
kax = kax(:);

% For the case where the root hydraulic conductivity is allowed to increase
% with depth, a linear increasing effect is considered
if rhc == 1     % if conductivity is allowed to increases with depth
    krad(rinds) = krad(rinds) + (z(rinds)/1000).*krad(rinds);   % linearly increasing 
    kax(rinds)  = kax(rinds) + (z(rinds)/1000).*kax(rinds);
end



% Root Pressure Potential --> HR Model
    % For the top soil layer
      j = 1;
      den    = z(j+1) - z(j);
      amx(j) = 0;
      bmx(j) = kax(j)/den + krad(j);
      cmx(j) = -kax(j)/den;
      dmx(j) = krad(j)*smp(j) - TR - kax(j);

    % For the middle soil layers
      for j = 2:N - 1
          den1   = z(j) - z(j-1);
          den2   = z(j+1) - z(j);
          amx(j) = -kax(j-1)/den1;
          bmx(j) = kax(j-1)/den1 + kax(j)/den2 + krad(j);
          cmx(j) = -kax(j)/den2;
          dmx(j) = krad(j)*smp(j) + kax(j-1) - kax(j);
      end 

    % For the bottom soil layer
      j = N;
      den    = z(j) - z(j-1);
      amx(j) = -kax(j-1)/den;
      bmx(j) = kax(j-1)/den + krad(j);
      cmx(j) = 0;
      dmx(j) = krad(j)*smp(j) + kax(j-1);


    % Solve for root pressure potential using tridiagonal matric solver
    aamx = amx(ginds);
    bbmx = bmx(ginds);
    ccmx = cmx(ginds);
    ddmx = dmx(ginds);
    nl = length(ginds);
%    rpp_ginds = TRIDIAG(nl, aamx, bbmx, ccmx, ddmx); 
    rpp = TRIDIAG(N, amx, bmx, cmx, dmx); 
    

% Weighted mean smp over root uptake profile [mm]
rpp(nrinds) = 0;
rpp = rpp(:);
rpp_weight = sum(rpp(rinds).*rootfr(rinds)/sum(rootfr(rinds)));% * mmH2OtoMPA;

