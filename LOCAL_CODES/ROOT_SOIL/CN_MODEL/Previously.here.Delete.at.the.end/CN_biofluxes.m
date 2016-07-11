% This function computes the fluxes of bioturbatino in each
% Horizon

function [Cin_m3, Cout_m3, Cin_m2, Cout_m2] = CN_biofluxes (Clsim, deltaz_diff, dz_diff, top, bottom,D)

nlayer = length(Clsim);

% allocate vectors
Cin_m2 = nan(nlayer,1);
Cout_m2 = nan(nlayer,1);

Cin_m2(1) =  top(1);
Cout_m2(1) = -(Clsim(2)-Clsim(1))*D/deltaz_diff(1);

for ii=2:1:nlayer-1
    Cin_m2(ii) = -(Clsim(ii)-Clsim(ii-1))*D/deltaz_diff(ii-1)+top(ii);
    Cout_m2(ii) = -(Clsim(ii+1)-Clsim(ii))*D/deltaz_diff(ii);
end
Cin_m2(nlayer) = -(Clsim(nlayer)-Clsim(nlayer-1))*D/deltaz_diff(nlayer-1);
Cout_m2(nlayer) = bottom;

Cin_m3 = Cin_m2./dz_diff;

Cout_m3 = Cout_m2./dz_diff;

