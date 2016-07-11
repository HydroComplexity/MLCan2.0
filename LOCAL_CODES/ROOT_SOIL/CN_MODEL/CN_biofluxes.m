% This function computes the fluxes of bioturbatino in each
% Horizon

function [Cin_m2, Cout_m2, difbio_m2, Bioflux] = CN_biofluxes (Clnew, dz, deltaz, D, BC)

nlayer = length(Clnew);

% allocate vectors
Bioflux = zeros(nlayer,2);
Cin_m2_m = zeros(nlayer,2);
Cout_m2_m = zeros(nlayer,2);
Cin_m2 = zeros(nlayer,1);
Cout_m2 = zeros(nlayer,1);

%*************************************************************************
% Generate matrix Bioflux with size [nsoil x 2], first column is the flux
% with information of flux from above and second column shows te flux 
% from below

% Top
D2 = D(1);
delz2 = deltaz(1);
Bioflux(1,1) =  BC(1);
Bioflux(1,2) = (Clnew(2)-Clnew(1))*D2/delz2;

% middle
for ii=2:1:nlayer-1
    D1 = D(ii-1);
    D2 = D(ii);
    delz1 = deltaz(ii-1);
    delz2 = deltaz(ii);
    Bioflux(ii,1) = -(Clnew(ii)-Clnew(ii-1))*D1/delz1;
    Bioflux(ii,2) = (Clnew(ii+1)-Clnew(ii))*D2/delz2;
end

% bottom
D1 = D(nlayer-1);
delz1 = deltaz(nlayer-1);
Bioflux(nlayer,1) =  -(Clnew(nlayer)-Clnew(nlayer-1))*D1/delz1;
Bioflux(nlayer,2) = -BC(nlayer);

%*************************************************************************
% Now compute the net inflow of carbon into layer 
difbio_m2 = sum(Bioflux,2);
Cin_m2_m(Bioflux>0) = Bioflux(Bioflux>0);
Cout_m2_m(Bioflux<0) = -Bioflux(Bioflux<0);

Cin_m2 = sum(Cin_m2_m,2);
Cout_m2 = sum(Cout_m2_m,2);


% % middle
% for ii=2:1:nlayer-1
%     D1 = D(ii-1);
%     D2 = D(ii);
%     delz1 = deltaz(ii-1);
%     delz2 = deltaz(ii);
%     Cin_m2(ii) = -(Clnew(ii)-Clnew(ii-1))*D1/delz1;
%     Cout_m2(ii) = -(Clnew(ii+1)-Clnew(ii))*D2/delz2;
% end
% 
% % bottom
% D1 = D(nlayer-1);
% delz1 = deltaz(nlayer-1);
% Cin_m2(nlayer) =  -(Clnew(nlayer)-Clnew(nlayer-1))*D1/delz1;
% Cout_m2(nlayer) = BC(nlayer);
% 
% 
% % convert units from m3 to m2
% Cin_m3 = Cin_m2./dz;
% Cout_m3 = Cout_m2./dz;
% 
