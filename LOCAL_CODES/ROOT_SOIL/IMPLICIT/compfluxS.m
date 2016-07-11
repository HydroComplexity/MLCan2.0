function [flux_s,flux_sr,flux_sr_all,ft,fb,dwat,flux_Ss, mberror,mberrormm]=...
compfluxS(psif,psiroot,kf,krad,ft,fb,dzsoi,zsoi,nl_soil, nstype,thetai,thetaf,dtime,nspecies,SS1,psii,thetaant)

%=========================================================================
% This function computes the water fluxes between layers in the soil colum
% and also it computes the fluxes between the soil and the root.   
%  
%
% Written by Juan Camilo Quijano, UIUC, 2008
%
%------------------------- Input Variables -------------------------------
%       psif            % [mm] Soil matrix solution
%       psiroot         % [mm] Root Water Potential
%       kf              % [mm/s] Hydraulic conductivity at interface between soil layers
%       krad            % Root Radial conductivities [1/s] 
%       ft              % [mm/s] Flux Top boundary condition.
%       fb              % [mm/s] Flux boundary condition at the botoom. 
%       dzsoi           % [mm] Vector with Layer thicknesses
%       zsoi            % [mm] Vector with depth of the different nodes in the layers  
%       nl_soil         % [] Number of layes
%       type            % 1.OK, 2. Supersaturation First Layer, 3. Super Saturation in between layers, thus hor_drainage >0 
%       thetai          % [] Volumetric water content in previous time step 
%       thetaf          % [] Solution of Volumetric Water Content in current time step
%       dtime           % [s] Time step
%       nspecies        % nspecies
%------------------------- Output Variables ------------------------------
%       flux_s          % [mm/s] Fluxes between soil layers
%       flux_sr         % [mm/s] Fluxes between soil and roots at all layers
%       flux_sr_all     % [mm/s] Fluxes between soil and roots at all layers
%       ft              % [mm/s] Flux Top boundary condition.
%       fb              % [mm/s] Flux boundary condition at the botoom.
%       dwat            % [] Change in volumetric water content
%       mberror         % [mm/s] Mass Balance error in Soil-Root solution.
%       mberrormm       % [mm] Mass Balance error in Soil-Root solution.
%-------------------------------------------------------------------------

[den]=createden(nl_soil,zsoi);
flux_s = nan(nl_soil,1);
flux_sr_all =nan(nl_soil,1);

%Ss
%flux_Ss = (Ss.*thetaf)./(thetas*dtime).*(psif-psii);    
flux_Ss = SS1.*thetaant.*(psif-psii);    

% fluxes in soil
for i=2:1:nl_soil
    flux_s(i)=-((psif(i)-psif(i-1))/(den(i-1)))*kf(i-1)+kf(i-1);    
    flux_s(nl_soil+1)=fb;
end
%for i=1:1:nl_soil
%    flux_sr(i)=(psif(i)-psiroot1(i))*krad1(i) + (psif(i)-psiroot2(i))*krad2(i);  
%end
flux_sr=zeros(nl_soil,1);
for i=1:1:nspecies
    flux_sr=(psif-psiroot(:,i)).*krad(:,i)+flux_sr;    
    flux_sr_all(:,i)=(psif-psiroot(:,i)).*krad(:,i);
end

if nstype(2)==1
    flux_s(1)=flux_s(2)+flux_sr(1)+(((thetaf(1)-thetai(1))*dzsoi(1))/dtime);
    ft=flux_s(1);
else
    flux_s(1)=ft;
end

flux_s=flux_s(:);
flux_sr=flux_sr(:);

dwat=thetaf-thetai; 

mberror=ft-sum(flux_sr)-fb-sum((thetaf-thetai).*dzsoi')/dtime - sum(flux_Ss.*dzsoi'); %Mass balance error in [mm/s]
mberrormm=mberror*dtime;
