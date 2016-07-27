clc
clear all
close all

% Hourly weather to half-hourly
%NumStep=17568;
NumStep=17520;
WEATHER=importdata('weather_2011_AndyMetdata.txt');

PPT_IN=WEATHER.data(:,5);
PPT_OUT=zeros(NumStep,1); % [mm]

Rg_IN=WEATHER.data(:,6);
Rg_OUT=zeros(NumStep,1); 

Ta_IN=WEATHER.data(:,7);
Ta_OUT=zeros(NumStep,1); 

RH_IN=WEATHER.data(:,8);
RH_OUT=zeros(NumStep,1); 

U_IN=WEATHER.data(:,9);
U_OUT=zeros(NumStep,1); 


% % Hourly data to half-Hourly
% year_blo=zeros(size(Si_Y,1)*2,1); % year
% Ta_blo=zeros(size(Si_Y,1)*2,1); % ['C]
% Pa_blo=zeros(size(Si_Y,1)*2,1); % [kPa]
% Rg_blo=zeros(size(Si_Y,1)*2,1); % [W/m2]
% PPT_blo=zeros(size(Si_Y,1)*2,1); % [mm]
% RH_blo=zeros(size(Si_Y,1)*2,1); % [mm]
% VPD_blo=zeros(size(Si_Y,1)*2,1); % [mm]
% U_blo=zeros(size(Si_Y,1)*2,1); % [mm]

for i=1:NumStep/2    
    PPT_OUT(i*2-1)=PPT_IN(i)/2; % [mm]
    PPT_OUT(i*2)=PPT_IN(i)/2;
    
    Rg_OUT(i*2-1)=Rg_IN(i); % [W/m2]
    Rg_OUT(i*2)=Rg_IN(i);
 
    Ta_OUT(i*2-1)=Ta_IN(i);
    Ta_OUT(i*2)=Ta_IN(i);

    RH_OUT(i*2-1)=RH_IN(i); % [%]
    RH_OUT(i*2)=RH_IN(i);

    U_OUT(i*2-1)=U_IN(i); % [m/s]
    U_OUT(i*2)=U_IN(i);
end

% http://cronklab.wikidot.com/calculation-of-vapour-pressure-deficit
SVP= 610.7.*10.^(7.5.*Ta_OUT./(237.3+Ta_OUT)); %[Pa]
VPD = ((100 - RH_OUT)./100).*SVP .*(0.001/1); % [Pa] to [kPa]

% http://keisan.casio.com/exec/system/1224579725
P0=101.325; % [kPa]
h=220;
T=Ta_OUT;
P=P0.*(1-(0.0065.*h)./(T+0.0065.*h+273.15)).^(5.257);