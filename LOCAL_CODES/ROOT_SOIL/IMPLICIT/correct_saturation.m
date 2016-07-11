function [flux_s,flux_sr,flux_sr_all,ft,fb,dwat,mberror,mberrormm, hor_drainage,hor_drainage_lay, psicom]=...
correct_saturation(psicom,psiroot,kf,krad,ft,fb,dzsoi,zsoi,nl_soil,...
type,thetai,thetaf,dtime,nspecies, pentry, thetas)

%=========================================================================
% This function computes the water fluxes between layers in the soil colum
% and also it computes the fluxes between the soil and the root.   
%  
%
% Written by Juan Camilo Quijano, UIUC, 2008
%
%------------------------- Input Variables -------------------------------
%       psicom          % Head Solution in each layer.
%       psiroot         % Vector with depth of midpoint of each layer [mm]
%       kf              % Vector with Layer thicknesses information [mm]
%       krad            % Vector with depth of boundaries of layers [mm] 
%       ft              % Number of layespsic
%       fb              % Flux top boundary condition 
%       dzsoi           % Vector with Layer thicknesses information [mm]
%       zsoi            % [mm] Vector with depth of the different nodes in the layers  
%       nl_soil         % Number of layes
%       type            % Type of top boundary condition.
%                         Type = 1, flux boundary condition.
%                         Type = 2, head boundary condition, first layer saturated .
%       thetai          % Soil moisture in previous time step
%       thetaf          % Solution of soil moisture in current time step
%       dtime           % Time step
%       nspecies        % Number of species
%       pentry          % [mm] minimum soil suction, i.e., soil potential at saturation 
%       thetas          % [] Soil moisture at saturation 
%
%------------------------- Output Variables ------------------------------
%       flux_s          % [mm/s] Fluxes between soil layers
%       flux_sr         % [mm/s] Fluxes between soil and roots at all layers
%       flux_sr_all     % [mm/s] Fluxes between soil and roots at all layers
%       ft              % [mm/s] Flux Top boundary condition.
%       fb              % [mm/s] Flux boundary condition at the botoom.
%       dwat            % [] Change in volumetric water content
%       mberror         % [mm/s] Mass Balance error in Soil-Root solution.
%       mberrormm       % [mm] Mass Balance error in Soil-Root solution.
%       hor_drainage    % [mm] Horizontal drainage.
%       psicom          % [mm] New Head Solution in each layer.
%
%-------------------------------------------------------------------------
dzsoi = dzsoi(:);
% Indicator of saturated layers
ind = psicom > pentry;
psicom(ind) = pentry(ind);

[den]=createden(nl_soil,zsoi);
flux_s = nan(nl_soil,1);
flux_sr_all =nan(nl_soil,1);
dwat =nan(nl_soil,1);


for i=2:1:nl_soil
    flux_s(i)=-((psicom(i)-psicom(i-1))/(den(i-1)))*kf(i-1)+kf(i-1);    
    flux_s(nl_soil+1)=fb;
end
%for i=1:1:nl_soil
%    flux_sr(i)=(psicom(i)-psiroot1(i))*krad1(i) + (psicom(i)-psiroot2(i))*krad2(i);  
%end
flux_sr=zeros(nl_soil,1);
for i=1:1:nspecies
    flux_sr=(psicom-psiroot(:,i)).*krad(:,i)+flux_sr;    
    flux_sr_all(:,i)=(psicom-psiroot(:,i)).*krad(:,i);
end

flux_s(1)=ft;

flux_s=flux_s(:);
flux_sr=flux_sr(:);

% Compute change in theta
for ii=1:nl_soil
   dwat(ii) = (flux_s(ii)-flux_s(ii+1)-flux_sr(ii))/dzsoi(ii)*dtime;    
end
dwat = dwat(:);
thetaf = thetai + dwat;
% Indicator of saturation in theta
indsat = thetaf > thetas;

% Compute horizontal drainage_should be the same as mberror
% Fixing the constant soil layer problem
%hor_drainage_lay = zeros(12,1);
hor_drainage_lay = zeros(nl_soil,1);
hor_drainage = sum((thetaf(indsat)-thetas(indsat)).*dzsoi(indsat))/dtime;
hor_drainage_lay(indsat) = ((thetaf(indsat)-thetas(indsat)).*dzsoi(indsat))/dtime;

% Update thetaf to be maximum to posority
thetaf(indsat) = thetas(indsat);
dwat = thetaf - thetai;


mberror=ft-sum(flux_sr)-fb-sum((thetaf-thetai).*dzsoi)/dtime; %Mass balance error in [mm/s]
mberrormm=(mberror - hor_drainage)*dtime;
