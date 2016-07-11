function [CNa, ADD, CNo, OUT, ADD_bio, ADD_ex, ADD_net] = CN_addlitter (FORCING, PARAMS,...
    VERTSTRUC, VARIABLES, CONSTANTS, SWITCHES, rootfr) 

% DE REFERENCE BLOCKS
% Dongkook Woo - Edit
if SWITCHES.CN.NupRootBiomass == 1
    TBMla = FORCING.TBMla; % Total input of litter in current time step [gr / m^2]
    CR_ratio = PARAMS.CN.CR_ratio;% factor that determines how much fraction of litter
                                  % goes below ground compare with above ground
    CNveg = PARAMS.CN.CNveg;
elseif SWITCHES.CN.NupRootBiomass == 0
    TBMla   = FORCING.BMla;      % Input of litter in current time step for each species [gr/m2]
    TBMrb   = FORCING.BMrb;
    CNabove = FORCING.CNabove;
    CNbelow = FORCING.CNbelow;
end
% Dongkook Woo - Edit End

if SWITCHES.CN.Bioturbation;   % include litter
    nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
    dz_mm = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm];%       dz_mm = grid depth [mm]
else
    nl_soil = PARAMS.nl_soil;%       nl_soil = # soil layers        
    dz_mm = VERTSTRUC.dzsmm;%       dz_mm = grid depth [mm]    
end
nspecies = PARAMS.CanStruc.nspecies; % number of species
timestep = VARIABLES.timestep; %timestep = Current time step
dtime = CONSTANTS.dtime;%       dtime = Time step [1800 s]
CNl_ini = VARIABLES.CNl;%       Initial CN in litter 

%DECLARE MATRICES 
ADD_r=zeros(nl_soil,nspecies);      % Addition of litter from the below-ground
ADD_l=zeros(nl_soil,nspecies);      % Addition of litter from the abovw-ground
ADDtotal = zeros(nl_soil,1);
CNadd = nan(1,nspecies);

%*************************************************************************
%                   COMPUTE THE LITTER TO ADD
%*************************************************************************
% IMPORRRRRTAAAAAANTTTT: ALL THE ADD FLUXES ARE COPUTED IN gr/m2/day 

% Above ground
% Dongkook Woo - Edit
if SWITCHES.CN.NupRootBiomass == 1
    ADDa = TBMla(:,timestep);                 % Total biomass Drop in this time step as litter from above [gr / m^2]
    ADDa = ADDa * 86400/dtime;                % Convert from [gr/m^2/dtime s] to [gr/m^2/d]
    % Below ground
    ADDb = ADDa.*CR_ratio';                     % Total biomass Drop in this time from below [gr / m^2/d]
elseif SWITCHES.CN.NupRootBiomass == 0
    ADDa = TBMla(:,timestep);                 % Total biomass Drop in this time step as litter from above [gr / m^2]
    ADDa = ADDa * 86400/dtime;                % Convert from [gr/m^2/dtime s] to [gr/m^2/d]
    % Below ground
    ADDb = TBMrb(:,timestep);                 % Total biomass Drop in this time step as litter from above [gr / m^2]
    ADDb = ADDb * 86400/dtime;                % Convert from [gr/m^2/dtime s] to [gr/m^2/d]
end
% Dongkook Woo - Edit End

for jj=1:1:nspecies
    % Dongkook Woo - Edit
    if SWITCHES.CN.NupRootBiomass == 1
        CNadd(jj) = CNveg{jj}(timestep);
    elseif SWITCHES.CN.NupRootBiomass == 0
        CNadd(jj) = CNabove(jj,timestep)'; 
        CNbadd(jj) = CNbelow(jj,timestep)';
    end
    % Dongkook Woo - Edit End
    % add litter from below ground
    ADD_r(:,jj)=ADDb(jj).*rootfr(:,jj)./(dz_mm/1000); 
    % add litter above ground is distributed uniformly in first layers     
    ADD_l(:,jj)=zeros(nl_soil,1); 
    ADD_l(1,jj)= ADDa(jj)/(dz_mm(1)/1000);
    % Supersposition of litter from both
    ADDtotal(:,jj)=ADD_r(:,jj)+ADD_l(:,jj);
end

% Remember all these fluxes are computed [gr/m3/day]. Daily time scale
ADDr = sum(ADD_r,2);  % compute total addition of litter [gr/m3/d]
ADDl = sum(ADD_l,2);  % compute total addition of litter [gr/m3/d]
ADD = ADDr + ADDl; % compute total addition of litter [gr/m3/d]
matCN = repmat(CNadd,nl_soil,1);  % compute matriz of CN from CNada

% compute new CN as (C1+C2/(C1/CN1 + C2/CN2))
CNa = ADD./sum(ADDtotal./matCN,2); 
% For the case when CNa is nan, does not matter because ADD
% is zero there is not flux. We just put a simple average
indnan = isnan(CNa);
CNa(indnan) = mean(CNadd);

if SWITCHES.CN.Bioturbation
    ADD_ex = ADDr;        % Addition of litter from external source  [gr/m3/d]
    % Dongkook Woo - Edit
    if SWITCHES.CN.NupRootBiomass == 1
        CNa_ex = CNa;        % CN in the litter added from external source
    elseif SWITCHES.CN.NupRootBiomass == 0
        CNa_ex = sum(CNbadd'.*ADDb/sum(ADDb));
        if isnan(CNa_ex)
            CNa_ex = mean(CNbadd);
        end
    end
    % Dongkook Woo - Edit End

    ADD_bio = VARIABLES.Cl_bio_in;   % Addition of litter from bioturbation source  [gr/m3/d]
    CNa_bio = VARIABLES.CN_bio_in;   % CN ratio in the bioturbation source
    OUT_bio = VARIABLES.Cl_bio_out;  % Extraction of litter from bioturbation source   [gr/m3/d]
    CNo_bio = VARIABLES.CN_bio_out;        % CN ration in the layer where is the out in bioturbation
    
    % Compute first the addition
    ADD = (ADD_ex + ADD_bio);       % Total ADD from litter and bioturbation
    CNa = (ADD)./((ADD_ex./CNa_ex) + (ADD_bio./CNa_bio));  % Weighted average of the CN ration in the ADD organic matter
    % Again Check for zeros. If ADD is equal to zero CNa is the same as before
    indadd = ADD>0;
    CNa(~indadd) = CNl_ini(~indadd);         

    % Compute second the extraction    
    OUT = OUT_bio;              % Total out from layer
    CNo = CNo_bio;              % CN ratio in the layer of extraction
    
    % Compute net
    ADD_net = ADD - OUT;        % [gr/m3/d]    
else
    CNo = nan(nl_soil,1);
    OUT = nan(nl_soil,1);
    ADD_ex = nan(nl_soil,1);    % nan [gr/m3/d]
    ADD_bio = nan(nl_soil,1);   % nan [gr/m3/d]
    ADD_net = nan(nl_soil,1);   % nan [gr/m3/d]
end       
