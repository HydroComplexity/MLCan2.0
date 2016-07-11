function [ rootfr roottr nl_root ] = ROOTDIST_LOGISTIC(zn, zh, dz, PARAMS, nspecies)
%
%   Calculate the logistic root fraction distribution
%
%   INPTUTS:
%       zn = node depths [m]
%       dz = layer thicknesses [m]
%       Z50 = depth to which 50% of root biomass exists
%       Z95 = depth to which 95% of root biomass exists
%
%   Written By: Darren Drewry

% De-reference Structure Values

for ii=1:1:nspecies
    maxrootdepth = PARAMS.Soil.maxrootdepth(ii);
    z50 = PARAMS.Soil.z50f(ii);
    z95 = PARAMS.Soil.z95f(ii);
    cc = 1.27875 / (log10(z50) - log10(z95));
    fr = -dz*cc/z50 .* (zn/z50).^(cc-1) .* (1+(zn/z50).^cc).^-2;
    % Cut off roots at maximum depth
    mind = find(zh>maxrootdepth, 1, 'first');
    fr(mind+1:end) = 0;
    
    fr = fr / sum(fr);

    z50 = PARAMS.Soil.z50t(ii);
    z95 = PARAMS.Soil.z95t(ii);
    cc = 1.27875 / (log10(z50) - log10(z95));    
    tr = -dz*cc/z50 .* (zn/z50).^(cc-1) .* (1+(zn/z50).^cc).^-2;
    % Cut off roots at maximum depth
    mind = find(zh>maxrootdepth, 1, 'first');
    tr(mind+1:end) = 0;
    tr = tr / sum(tr);
    
    if isempty(mind)
        nl_root(ii)=length(fr);
    else
        nl_root(ii)=mind; 
    end
    rootfr(:,ii)=fr;
    roottr(:,ii)=tr;
end
