
function [Cl, VARIABLES] = CN_bioturbation (PARAMS, VARIABLES, CONSTANTS, FORCING, VERTSTRUC, SWITCHES, Cl, fTd)

% ALLOCATE MATRICES TO USE
nspecies = PARAMS.CanStruc.nspecies;
CNveg = nan(1,nspecies);
%   INPUTS:
% DE REFERENCE BLOCKS
%  VARIABLES structure
timestep = VARIABLES.timestep;          % timestep = Current time step
DECl = VARIABLES.DECl;
BD = VARIABLES.BD;
litter_dz = VARIABLES.SOIL.litterthickness;%   layer thickness [m]
CNl = VARIABLES.CNl;
%  FORCING structure
TBMla = FORCING.TBMla(:,timestep);      % Input of litter in current time step for each species [gr/m2]
ADDa = sum(TBMla);                      % Total input of litter in current time step [gr / m^2]  
%  PARAMS structure
for ii=1:1:nspecies
   CNveg(ii) = PARAMS.CN.CNveg{ii}(timestep); 
end
CNveg_mean = sum(CNveg'.*TBMla/sum(TBMla));
kbio = PARAMS.kbio;
nl_soil=PARAMS.nl_soil;                 % Number of layers in the soil
Clitter = PARAMS.CN.Clitter;
layerbio = PARAMS.CN.layerbio;
D= PARAMS.CN.D;
fHorO = PARAMS.CN.fHorO;
fTdl = PARAMS.CN.fTdl(timestep);
%  VERTSTRUC structure
dz = VERTSTRUC.dzs;
nz = VERTSTRUC.zns;                     % nodes depth in the grid [m]                                        
%  CONSTANTS structure
dtime = CONSTANTS.dtime;                % dtime = Time step [1800 s]

% BIOTURBATION AND LITTER FRAGMENTATION FLUX FROM LITTER TO O HORIZON 
%**************************************************************************
%PotBiot = PotBiot/(365*3600*24)*dtime;   %[grC/m2] in the time step 
Biotflux = Clitter(1)*litter_dz*kbio*fTdl; %[grC/m2/dtime] in the time step 
Cflux = min(Biotflux,Cl(1)*(litter_dz-0.005));     %[grC/m2/dtime] in the time step 
VARIABLES.biorate = Cflux;                % Bioturbation and fragmentation rate from litter
                                          % to horizon O. [gr/m2/dtime]

% change in total carbon in litter per m2
decom = (DECl(1)-BD(1))*litter_dz*dtime/86400;  % [grC/m2/dtime]
%decom = (Clitter-Cl(1))*litter_dz;    % just to check both should be the same

dClit_m2 = ADDa - Cflux - decom;    %[grC/m2]

% Conmpute new CN ratio in litter layer
CNl(1) = (Cl(1)*litter_dz*CNl(1) + dClit_m2*CNveg_mean)/(Cl(1)*litter_dz + dClit_m2);

% change in litter thickness
litterthickness = (Cl(1)*litter_dz + dClit_m2)/Clitter;  % [m] Satisfy constant litter concentration
Cl(1) = Clitter;                                          % Set litter concentration
deltathick = litterthickness - VARIABLES.SOIL.litterthickness; % [m]
VARIABLES.SOIL.litterthickness = litterthickness;    % [m]


%               BIOTURBATION DUE TO DIFFUSION 
%**************************************************************************
% Organize the layers that are needed for the calculation of diffusion
% Diffusion is computed until a given layer 'layerbio'
%ind_join = [1 1 2 2 3 4 5 6 6 7 7 7]'; 
ind_join = [1 2 3 4 5 6 7 8 9 10 11 12]'; 
[Cl_diff, dz_diff, zh_diff, zn_diff, deltaz_diff] = CN_layersdif (dz, Cl(2:nl_soil+1), layerbio, 1, ind_join); % for Carbon
[CNl_diff, dz_diff, zh_diff, zn_diff, deltaz_diff] = CN_layersdif (dz, CNl(2:nl_soil+1), layerbio, 1, ind_join); % for CNratio


Bio_Top = zeros(length(Cl_diff),1);
Bio_Top(1) = Cflux*(1-fHorO);       % Bioturbation rate 1st layer [grC/m2/dtime]
Bio_Top(2) = Cflux*fHorO;       % Bioturbation rate 2nd layer [grC/m2/dtime]

% change units in bioturbation difussion constant
D = D/100/100/365/48;          % from [cm2/year] to [m2/dtime]

%Compute the matrices
ndim = length(Cl_diff);
[A, B]  = CN_biomatrices  (dz_diff, deltaz_diff, D, 1, Bio_Top, ndim);

% Compute new states of C
Clsim = A^(-1)*(Cl_diff+B);
Clsim = Clsim(:);
% Compute fluxes in and fluxes out in each layer
[Cin_m3, Cout_m3, Cin_m2, Cout_m2] = CN_biofluxes (Clsim, deltaz_diff, dz_diff, Bio_Top, 0, D);
 
difbio = Clsim-Cl_diff;     % [gr/m3/dtime] 

%**************************************************************************

%       ALLOCATE THE DIFFERENCE IN CARBON IN EACH LAYER
dzinput(:,1) = dz(1:layerbio);
dzinput(1:length(dz_diff),2) = dz_diff;

% Organize the ouput in layers format using the information resulting from
% bioturbation which is in the horizons format
%ind_join = [1 1 2 2 3 4 5]'; 
ind_join = [1 2 3 4 5 6 7]'; 
[Cl_change, dz_dif, zh_dif, zn_dif, deltaz_dif] = CN_layersdif (dzinput, difbio, layerbio, 2, ind_join);
[Cl_in, dz_dif, zh_dif, zn_dif, deltaz_dif] = CN_layersdif (dzinput, Cin_m3, layerbio, 2, ind_join);
[Cl_out, dz_dif, zh_dif, zn_dif, deltaz_dif] = CN_layersdif (dzinput, Cout_m3, layerbio, 2, ind_join);

% Compute the final output considering the litter layer and the layers that
% are not affected by bioturbation
Cl_bio_change = [0 ; Cl_change ; zeros(nl_soil - layerbio,1)];    % [gr/m3/dtime] 
Cl_bio_in = [0 ; Cl_in ; zeros(nl_soil - layerbio,1)];    % [gr/m3/dtime]
Cl_bio_out = [0 ; Cl_out ; zeros(nl_soil - layerbio,1)];    % [gr/m3/dtime]

VARIABLES.Cl_bio_change = Cl_bio_change*86400/dtime;                % [gr/m3/d]
VARIABLES.Cl_bio_in = Cl_bio_in*86400/dtime;                        % [gr/m3/d]                    
VARIABLES.Cl_bio_out = Cl_bio_out*86400/dtime;                      % [gr/m3/d]

% Organize vectors of CNl.
% Compute the vector of CN for Input and Output in each layer that is affected by 
% bioturbation to calculate the new CN later in other function
CN_bio_in = [CNveg_mean ; CNl(1) ; CNl(1) ; CNl_diff(1) ; CNl_diff(1) ; CNl_diff(2) ; CNl(6:nl_soil)];
VARIABLES.CN_bio_in = CN_bio_in;
CN_bio_out = [CNl(1) ; CNl_diff(1) ; CNl_diff(1) ; CNl_diff(2) ; CNl_diff(2); CNl(5:nl_soil)];
VARIABLES.CN_bio_out = CN_bio_out;

%**************************************************************************
% CHECK MASS BALANCE BY SOLUTION OF BIOTURBATION

% CARBON
input = Cflux;     % [grN/m2/dtime]
Cchange = (Cl_bio_in - Cl_bio_out).*[litterthickness ; dz]; %[grN/m2/dtime]
VARIABLES.bioCerror = (input-sum(Cchange))*86400/dtime;  %[grC/m2/d] 


% NITROGEN
input = Cflux/(CNl(1));     % [grN/m2/dtime]
Nchange = (Cl_bio_in./CN_bio_in - Cl_bio_out./CN_bio_out).*[litterthickness ; dz]; %[grN/m2/dtime]
VARIABLES.bioNerror = (input-sum(Nchange))*86400/dtime;  %[grN/m2/d] 


