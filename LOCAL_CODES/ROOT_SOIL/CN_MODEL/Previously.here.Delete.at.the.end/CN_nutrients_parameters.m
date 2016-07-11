
% De reference block
TBMla = FORCING.TBMla;
CNh = PARAMS.CN.CNh;
rr = PARAMS.CN.rr;
fc= 0.39;     %Saxton et al (1986)    
%************************************************************************
% INTERPOLATE SOIL MOISTURE AND SOIL TEMPERATURE

% SOIL MOISTURE
zblo = [0 10 30 50 1000]; % depth of measurements on blodgett 
% 0 and 5 cm
ind=~isnan(SWC1_in);
theta_blo(1) = mean(SWC1_in(~isnan(SWC1_in)))/100;
theta_blo(2) = mean(SWC1_in(~isnan(SWC1_in)))/100;
% 10 cm
ind=~isnan(SWC2_in);
theta_blo(3) = mean(SWC2_in(~isnan(SWC2_in)))/100;
% 30 cm  
ind=~isnan(SWC3_in);
theta_blo(4) = mean(SWC3_in(~isnan(SWC3_in)))/100;
% 1000 close to saturation field capacity
theta_blo(5) = 0.3;

% interpolate
thetaave = interp1(zblo, theta_blo, zhsmm/10);

% TEMPERATURE 

% SOIL MOISTURE
zblo = [0 5 7.5 10 30 50 1000]; % depth of measurements on blodgett 

% 0 and 5 cm
ind=~isnan(Ts1_in);
T5v(1) = mean(Ts1_in(ind));
ind=~isnan(Ts4_in);
T5v(2) = mean(Ts4_in(ind));
ind=~isnan(Ts7_in);
T5v(3) = mean(Ts7_in(ind));
temp_blo(1) = mean(T5v);
temp_blo(2) = mean(T5v);

% 7.5 cm
ind=~isnan(Ts2_in);
T7v(1) = mean(Ts2_in(ind));
ind=~isnan(Ts5_in);
T7v(2) = mean(Ts5_in(ind));
temp_blo(3) = mean(T7v);

% 10 cm
ind=~isnan(Ts3_in);
T10v(1) = mean(Ts3_in(ind));
ind=~isnan(Ts6_in);
T10v(2) = mean(Ts6_in(ind));
ind=~isnan(Ts8_in);
T10v(3) = mean(Ts8_in(ind));
temp_blo(4) = mean(T10v);

% 30 cm
ind=~isnan(Ts9_in);
T30v = mean(Ts9_in(ind));
temp_blo(5) = T30v;

% 50 and 1000 cm
ind=~isnan(Ts10_in);
T50v = mean(Ts10_in(ind));
temp_blo(6) = T50v;
temp_blo(7) = T50v;

tempave = interp1(zblo, temp_blo, zhsmm/10);

%*************************************************************************
% COMPUTE STEADY STATE CONCENTRATIONS FROM DATA

zblo = [0      2.5   10.5   23    40    65   155 400 1000]; % depth of measurements on blodgett 
Cl_blo =   [45000 45000 40000 30000  14400 7700 3906 474  10 ];
Cb_blo =   [738    738   396    129   30    10    4   1    1 ]; 

% interpolate
Cl_bloi = interp1(zblo, Cl_blo, zhsmm/10);
Cl_bloi = Cl_bloi(:);
Cb_bloi = interp1(zblo, Cb_blo, zhsmm/10);
Cb_bloi = Cb_bloi(:);
Ch_bloi = ones(12,1);

Cl_bloii=45000*exp(-0.023*(zhsmm/10-2.5));

Cl_bloii(1:8) = Cl_bloii(1:8); 
%************************************************************************
% COMPUTE FACTORS

% soil moisture (base on Porporato soil moisture)
rwc = (thetaave-theta_dry)./(porsl-theta_dry);
if rwc<fc
    fSd=rwc./fc;
else
    fSd=fc./rwc;
end    

% temperature 
indT = tempave > 25;
fTd(indT) = 1;              
fTd(~indT) = 1.9.^((tempave(~indT)-25)/10); 
fTd = fTd(:);

%************************************************************************
% COMPUTE ADD LITTER
%TBMla = infos.TBMla;
cumlit = cumsum(TBMla,2);
ADDLEAF = cumlit(:,end)/365;      % [gC/m^2/d]
ADDROOT = ADDLEAF*PARAMS.CN.CR_ratio;      % [gC/m^2/d]
ADDLEAFtem = repmat(ADDLEAF', 12,1);
ADDROOTtem = repmat(ADDROOT', 12,1);
dzsmmtem=repmat(dzsmm,1,nspecies);
ADD2v = ADDROOTtem.*rootfr./(dzsmmtem/1000);
ADD2 = sum(ADD2v,2);

ADD1v=zeros(size(ADD2v));
ADD1v(1,:) = ADDLEAF'./(dzsmmtem(1,:)./1000);

fHorO = dzsmm(2)/(dzsmm(2)+dzsmm(1)); %****** same concentratio to O horizon
if SWITCHES.CN.Bioturbation  % MEANS ALSO BIOTURBATION    
    ADD1v(2,:) = ADDLEAF'*fHorO./(dzsmmtem(2,:)./1000);
    ADD1v(1,:) = ADD1v(1,:) - ADDLEAF'*fHorO./(dzsmmtem(1,:)./1000);
end
ADD1 = sum(ADD1v,2);
ADD=ADD1+ADD2;
%Compute the Total Additio in the year and the Potential Bioturbation
totalyearADD = sum(sum(TBMla));    %[grC/m2/year]
PotBiot = totalyearADD * PARAMS.CN.bio.biofactor; %[grC/m2/year]
VARIABLES.CN.PotBiot = PotBiot;       %[grC/m2/year]
    

%************************************************************************
% COMPUTE CONSTANTS
if SWITCHES.comcontants

    phi=1; %Assume phi=1 to compute constants
    rhc = min (0.25);  
    
% Compute fluxes due to bioturbation
    layerbio = PARAMS.layerbio;
    D= PARAMS.CN.D;
    D = D/100/100/365;          % from [cm2/year] to [m2/day]
    if SWITCHES.CN.Bioturbation
        % Create vector that groups layers for bioturbation
        indvec = [1 1 2 3 4 5 6 7 8 9 10 11];
        %indvec = [1 2 3 4 5 6 7 8 9 10 11 12];
        % Create new vectors for dz and hz for new grouped layers
        dzmm_cum = accumarray(indvec', dzsmm);
        zhmm_cum = cumsum(dzmm_cum);
        % Based on this new vectors define length and distance btn nodes
        ncum = length(zhmm_cum);        
        deltaz = (zhmm_cum(2:ncum) - zhmm_cum(1:ncum-1))/1000; % Delta z btn nodes in [m]       
        % Compute concentration for new grouped layers using weighted
        % average base on each initial (sublayer) layer thickness
        Cl_bloi_m2 = Cl_bloi.*(dzsmm/1000);
        Cl_blo_m2_cum = accumarray(indvec', Cl_bloi_m2);
        Cl_blo_cum = Cl_blo_m2_cum./(dzmm_cum/1000);
        % Compute the fluxes due to Bioturbation using the grouped layers
        Fluxbiot_cum = -D*(Cl_blo_cum(2:ncum)-Cl_blo_cum(1:ncum-1))./deltaz; % Flux btn layers 
        Fluxbiot_cum = [Fluxbiot_cum ; 0];  % No flux in bottom layer        
        Fluxbiot_cum(indvec(layerbio):ncum) = 0;  % put zeros after cutt
        Fluxbiot_in_cum = [0 ; Fluxbiot_cum(1:ncum-1)];  %No flux in bottom layer        
        Fluxbiot_out_cum = [Fluxbiot_cum];
        deltabiot_cum = (Fluxbiot_in_cum - Fluxbiot_out_cum)./(dzmm_cum/1000);
        
        f_dist = dzsmm./dzmm_cum(indvec);
        deltabiot_cum_m2 = deltabiot_cum.*(dzmm_cum/1000);
        deltabiot_m2 = deltabiot_cum_m2(indvec).*f_dist;
        deltabiot = deltabiot_m2./(dzsmm/1000);
        
        ADD = ADD + deltabiot;
    end
        for ii=1:1:nl_soil
          [kd(ii), kl(ii), kh(ii)] = ...
            CN_COMCONSTANTS(ADD(ii), Cb_bloi(ii), Cl_bloi(ii), Ch_bloi(ii), fSd(ii), fTd(ii), rhc, rr, phi, SWITCHES);
        end
else
        load constants.mat;
end
    if SWITCHES.CN.Bioturbation   % if Bioturbation is on, litter layers is included
        kd=[kd(1); kd(:)];
        PARAMS.CN.kd = kd;
        kl=[kl(1); kl(:)];
        PARAMS.CN.kl = kl;
        kh=[kh(1); kh(:)];
        PARAMS.CN.kh = kh;
    else        
        kd=kd(:);
        PARAMS.CN.kd = kd;
        kl=kl(:);
        PARAMS.CN.kl = kl;
        kh=kh(:);
        PARAMS.CN.kh = kh;            
    end
    
