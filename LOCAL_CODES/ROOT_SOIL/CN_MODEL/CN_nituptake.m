function [CC_up_m3, CC_up_m2,CC_up_all_m2] = CN_nituptake (PARAMS, rootarq,...
    CCdd, nl_soil, dtime, ioncase)


% Dereference block
nspecies = PARAMS.CanStruc.nspecies; % number of species;
diam = PARAMS.Soil.rootarq.diam;       % [mm] 
if ioncase == 1 
    vmax = PARAMS.vmaxnit;
    km = PARAMS.kmnit;
    % change units
    vmax = vmax*10^(-6)*14*24;  % from [mumol_no3/(g_rr hr)] to [g_N/g_rr day] 
    km = km*10^(-3)*14; % from [mmol_no3/(m3)] to [g_N/m3]
elseif ioncase == 2
    vmax = PARAMS.vmaxamm;
    km = PARAMS.kmamm;
    % change units
    vmax = vmax*10^(-6)*14*24;  % from [mmol_NH4/(g_rr hr)] to [g_N/g_rr day] 
    km = km*10^(-3)*14; % from [mmol_NH4/(m3)] to [g_N/m3]    
end
RMId = rootarq.RMId;   % [g_root/m2]
RVID = sum(rootarq.RVID,2);  % [V_rZD/m2]

% Allocate memory
CC_up_all_m2 = zeros(nl_soil,nspecies);

for ii=1:1:nspecies
    flux = zeros(nl_soil,1);
    for jj=1:1:length(diam)
        flux = (vmax.*CCdd.*RMId(:,ii,jj))./(km+CCdd) + flux;  %[gr N/m2/d]
    end
    CC_up_all_m2(:,ii) = flux;   %[gr N/m2/d]
end
CC_up_m2 = sum(CC_up_all_m2,2); %[g/m2/d]
% Compute the fraction of CC taken btn species
frcspc = CC_up_all_m2./repmat(CC_up_m2,1,nspecies);
%**** correct nan in frcspc  ************
indnan = isnan(frcspc);
frcspc(indnan) = 0; 
%**************************************************************************

total_rr = CC_up_m2./RVID*dtime/86400; %[g/m3/d]
ind_mb = CCdd < total_rr;        % Check mass availability in that layer in [gr/m^3]
CC_up_m2(ind_mb) = CCdd(ind_mb).*RVID(ind_mb)*86400/dtime;    % [gr/m^2/d]
CC_up_m3 = CC_up_m2./RVID;     %[g/m3/d]

%**** correct nan in frradvectionCCdd by zeros in RVID ************
indnan = isnan(CC_up_m3);
CC_up_m3(indnan) = 0; 
