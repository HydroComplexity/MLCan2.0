function [rrdiff_m2, rrdiffCC_m3, rrdiffCCdd_m3] =...
    CN_diffwu (PARAMS, rootarq, nl_soil, dz, nz, zh, CC,...
    CCdd, sm, dtime, De)  

% De Diffusivity [m2/d]

% Dereference block
diam = PARAMS.Soil.rootarq.diam;  %[mm]
deltarZD = PARAMS.Soil.rootarq.deltarZD;  %[mm]
RAIDd=rootarq.RAId;     % [m2 leaf / m2 ground]
RVID = sum(rootarq.RVID,2);
%jaeger = PARAMS.CN.jaeger;  


% Allocate needed variables
flux = zeros (nl_soil,1);
out_rrdiff = zeros(nl_soil,1);
in_rrdiff = zeros(nl_soil,1);
maxflux = zeros(nl_soil,1);
rrdiff_m2 = zeros(nl_soil,1);

RAIdT = squeeze(sum(RAIDd,2)); 

for ii =1 :1:length(diam)
    ra = (diam(ii)/2 + deltarZD)/1000;          % radiuous of disturbance zone in m
%    ra = (diam(ii)/2)/1000;          % radiuous of disturbance zone in m    
%    val = De.*(dtime/86400)./(ra.^2);    % (D*dtime/(r^2)) Values to compute function in Jaeger
%    ff = interp1(jaeger(1:26,1),jaeger(1:26,2),val,'linear','extrap');  % compute function Jager and Clark (1942) 
%    flux = ((CC - CCdd).*De/ra).*RAIdT(:,ii).*ff + flux;  %[g/m2/d]
    flux = ((CC - CCdd).*De/ra).*RAIdT(:,ii) + flux;
end

% 1. Compute maximum flux can be acchieved until equilibrium in dtime
%**************************************************************************
wtind = CC > CCdd;
wtnind = ~wtind;
% VERRYYY IMPORTANT THE MAXIMUM FLOW IS THE DIFFERENCE IN CONCENTRATION
% TIME THE RVID FOR BOTH CASES BCAUSE ZD IS THE LIMITING SOURCE OR SINK
maxflux(wtind) = (CC(wtind)-CCdd(wtind)).*RVID(wtind)*86400/dtime;
maxflux(wtnind) = (CCdd(wtnind)-CC(wtnind)).*RVID(wtnind)*86400/dtime;

%[g/m2/d]  Again as in advection negative out positive into soil
rrdiff_m2(wtind) = -(min(abs(flux(wtind)),maxflux(wtind)));
rrdiff_m2(wtnind) = (min(abs(flux(wtnind)),maxflux(wtnind)));
%**************************************************************************
     
%2. check soil availability
%**************************************************************************
out_rrdiff(wtind) = rrdiff_m2(wtind)./dz(wtind)*dtime/86400;
out_rrdiff(~wtind) = 0;       
ind_mb = CC < abs(out_rrdiff);        % Check mass availability in that layer in [gr/m^3]
rrdiff_m2(ind_mb) = -CC(ind_mb).*dz(ind_mb)*86400/dtime;    % [gr/m^2/d]

% check soil availability
in_rrdiff(wtnind) = rrdiff_m2(wtnind)./(RVID(wtnind))*dtime/86400;
in_rrdiff(~wtnind) = 0;       
ind_mb = CCdd < abs(in_rrdiff);        % Check mass availability in that layer in [gr/m^3]
rrdiff_m2(ind_mb) = CCdd(ind_mb).*RVID(ind_mb)*86400/dtime;    % [gr/m^2/d]
%**************************************************************************

rrdiffCC_m3 = rrdiff_m2./dz;            % [gr/m^3/d]
rrdiffCCdd_m3 = -rrdiff_m2./RVID;            % [gr/m^3/d]
        
%**** correct nan in frradvectionCCdd by zeros in RVID ************
indnan = isnan(rrdiffCCdd_m3);
rrdiffCCdd_m3(indnan) = 0; 
     