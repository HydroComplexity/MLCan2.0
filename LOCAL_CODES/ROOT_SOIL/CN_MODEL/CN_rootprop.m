function [rootarq] = CN_rootprop (FORCING, PARAMS, VARIABLES, VERTSTRUC, SWITCHES)
% INPUT INFORMATION
%diam [mm] Column vector with information of mean diameter in each class
%SRL [m/g] Column vector with information of Specific root length for each class
%dist []  Column vector with distribution of length in each class. sum(dist) = 1
%deltar [mm] Distance that defines the zone of disturbance

% OUTPUT INFORMATION
%LRI [m/m2]. Length Root Index. Root length per ground area
%ARI [g/m2]. Area Root Index. Surface area of root per ground area
%VRI [g/m2]. Volumetric Root Index. Volume of root per ground area
%ARId [g/m2]. Area Root Index in the DZ. Total surface area of the disturbance zone per ground area
%VRI [g/m2]. Volumetric Root Index in the DZ. Volume of the DZ per ground area



%  DE REFERENCE BLOCKS
nspecies = PARAMS.CanStruc.nspecies; % number of species
rootfr = VERTSTRUC.rootfr;
if SWITCHES.CN.Bioturbation;
   nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers  
   % Dongkook Woo - Edit
   %rootfr  =[zeros(1:nspecies); rootfr];
   rootfr  =[zeros(1,nspecies); rootfr];
   % Dongkook Woo - Edit End
else
   nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
end

 
% Root arquitecture
    diam = PARAMS.Soil.rootarq.diam;       % [mm] 
    dist = PARAMS.Soil.rootarq.dist;       % 
    SRL = PARAMS.Soil.rootarq.SRL;        % [m/g]
    deltarZD = PARAMS.Soil.rootarq.deltarZD;       % [mm]
    tt = VARIABLES.timestep;               % [s] 
% convert diameters from [mm] to [m]
    diam = diam/1000;
    deltarZD = deltarZD/1000;

    
% ALLOCATE SPACES
RLId = zeros(nl_soil,nspecies,length(dist));
mass = zeros(1,nspecies);    
%mattem = nan(nl_soil, length(diam), nspecies);
    
for ii=1:nspecies
    mass(ii) = FORCING.LMI(ii,tt);   % [g]
end

% COMPUTATION OF RMI, RLI, RAI, RVI, RAID, RVID
RMI = repmat(mass,nl_soil,1).*rootfr;


% computation of RLI, RMI, RAI and RVI
% 1 RLI
RLI = RMI./(sum(dist./SRL));
for kk=1:1:length(dist)
    RLId(:,:,kk) = RLI*dist(kk);
end

% 2 RMI 
%compute first matrix (di/SRL))
vectem = dist./SRL;
mattem = repmat(vectem, [nl_soil,1,nspecies]);
RMId = permute(mattem, [1 3 2]);
%compute now L for each 
RLIsp = repmat(RLI,[1,1,4]);
RMId = RMId.*RLIsp;
%CHECK sum(RMId,3) = RMI


% 3. RAI 
%compute first matrix (2di(ri^2))
vectem = dist.*diam;
mattem = repmat(vectem, [nl_soil,1,nspecies]);
RAId = permute(mattem, [1 3 2]);
%compute now L for each 
RLIsp = repmat(RLI,[1,1,4]);
RAId = RAId.*RLIsp*pi;
RAI = RLI*pi*sum(dist.*diam);
%CHECK sum(RAId,3) = RAI

% 4. RVI
%compute first matrix (di(ri^2))
vectem = dist.*(diam/2).^2;
mattem = repmat(vectem, [nl_soil,1,nspecies]);
RVId = permute(mattem, [1 3 2]);
%compute now L for each 
%RLIsp = repmat(RLI,[1,1,4]);
RVId = RVId.*RLIsp*pi;
RVI = RLI*pi*sum(dist.*(diam/2).^2);
%CHECK sum(RVId,3) = ARI

% 5. RAI disturbance
%compute first matrix (2di(ri + deltar))
vectem = 2*dist.*(diam/2+deltarZD);
mattem = repmat(vectem, [nl_soil,1,nspecies]);
RAIDd = permute(mattem, [1 3 2]);
%compute now L for each 
%RLIsp = repmat(RLI,[1,1,4]);
RAIDd = RAIDd.*RLIsp*pi;
RAID = RLI*pi*sum(2*dist.*((diam/2)+deltarZD));
%CHECK sum(RAIDd,3) = RAID


% 6. RVI disturbance
%compute first matrix (di(2rideltar + deltar^2))
vectem = dist.*(2*(diam/2)*(deltarZD)+deltarZD^(2));
mattem = repmat(vectem, [nl_soil,1,nspecies]);
RVIDd = permute(mattem, [1 3 2]);
%compute now L for each 
%RLIsp = repmat(RLI,[1,1,4]);
RVIDd = RVIDd.*RLIsp*pi;
RVID = RLI*pi*sum(dist.*(2*(diam/2)*deltarZD+(deltarZD^2)));
%CHECK sum(RVIDd,3) = RVID


% STORE INFORMATION IN STRUCTURES

rootarq.RMI = RMI;
rootarq.RMId = RMId;
rootarq.RAId = RAId;
rootarq.RAI = RAI;
rootarq.RLId = RLId;
rootarq.RLI = RLI;
rootarq.RVId = RVId;
rootarq.RVI = RVI;
rootarq.RAIDd = RAIDd;
rootarq.RAID = RAID;
rootarq.RVIDd = RVIDd;
rootarq.RVID = RVID;

