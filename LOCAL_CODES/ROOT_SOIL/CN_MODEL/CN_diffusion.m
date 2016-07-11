function [DIFF_m3, DIFF_m2, TLCDFF_m2, TLCDFF_m3, De] = CN_diffusion(nl_soil, dz, nz, zh, CC, sm, dtime, ioncase)

%*************************************************************************
% Declare matrix and variables needed
Fo_diff = zeros (nl_soil,2);
DIFF_m2 = zeros (nl_soil,1);

soildif = 2;

%*************************************************************************
%       Compute Diffusivities
% Compute Effective Diffusivities based on soil moisture at each node
De = zeros(nl_soil,1);
Den = zeros(nl_soil-1,1);
switch soildif 
    case 1
        if ioncase == 1
            De = (0.36*+0.22.*sm)*10^(-6);  %[cm^2/s]
        else
            De = (0.98 + 11.06.*sm)*10^(-7); %[cm^2/s]
        end    
    case 2
        if ioncase == 1
            De = (1.41-19.9.*sm+114*(sm).^2)*10^(-6); %[cm^2/s]
        else
            De = (-0.87+12.2.*sm)*10^(-7); %[cm^2/s]
        end
end    
% check negative conductivities
indneg = De<0;
De(indneg) = 0;
% Compute the effectivie diffusivities at the interface layer using a 
% series approximation

% Dongkook Woo - Edit
%for ii=1:1:11
for ii=1:1:nl_soil-2
    % compute diffusivity in mieddle point btn layer ii and layer ii+1
    Den(ii) = (nz(ii+1)-nz(ii))/((zh(ii)-nz(ii))/De(ii) + (nz(ii+1)-zh(ii))/De(ii+1)) ; %[cm^2/s]
end
% Assume diffusivity in last interface is the same as diffusivity in bottom
% layer
%Den(12) = De(12);  %[cm^2/s]
Den(12) = De(nl_soil-1);  %[cm^2/s]
% Dongkook Woo - Edit End

%  Change units from [cm^2/s]  to [m^2/d]
Den = Den * 86400/((100)^2);
De = De * 86400/((100)^2); 
%*************************************************************************
%       Compute Fluxex


% Define Fo_adv which is a [nlayer X 2] matrix with out fluxes from each
% layer. 1st column is upward flow while 2nd column is downward

% Compute matrix from ii=1 up to ii=12. It computes the fluxes under each 
% layer ii
for ii = 1:1:nl_soil-1
    flux = ((CC(ii) - CC(ii+1))/(nz(ii+1)-nz(ii)))*Den(ii);  % [g/m2/d]
    if flux > 0
        Fo_diff(ii,2) = flux;         % [g/m2/d]
    else 
        Fo_diff(ii+1,1) = abs(flux);       % [g/m2/d]
    end
end    
    
% Assume no diffusion flux at the botton
Fo_diff(nl_soil,2) = 0;              % [g/m2/d]


% Compute total fluxex and check there is enough mass in each layer to
% satisfy that flux
F_diff = sum(Fo_diff,2);       % [gr/m^2/d]   
out_diff = F_diff./dz * dtime/86400;    % Total out due to advection in [gr/m^3/dtime]
ind_mb = CC < out_diff;        % Check mass availability in that layer in [gr/m^3]
F_diff(ind_mb) = 0;           % Those layers with no available are set to 0 flux  
Fo_diff(ind_mb,:) = 0;     
    
% Using the output fluxes matrix Fo_diff now we compute the flux into each layer 
DIFF_m2 = DIFF_m2 - F_diff;% Substract the fluxes out from each layer first   
% Add the fluxes from each layer from the 2 to the nlsoil -1
for ii = 2:nl_soil-1                
    DIFF_m2(ii-1) = DIFF_m2(ii-1) + Fo_diff(ii,1);  %[g/m^2/d]
    DIFF_m2(ii+1) = DIFF_m2(ii+1) + Fo_diff(ii,2);  %[g/m^2/d] 
end
% Now compute for the top  
DIFF_m2(2) = DIFF_m2(2) + Fo_diff(1,2);             %[g/m^2/d]
% and botom
DIFF_m2(nl_soil-1) = DIFF_m2(nl_soil-1) + Fo_diff(nl_soil,1); %[gr/m^2/d] 
TLCDFF_m2 = Fo_diff(nl_soil,2);           % Total leaching [gr/m^2/d]
%
%********************** WE CAN CHECK MASS BALANCE. COMPUTING THIS *****
%sum(DIFF_m2) + TLCDFF_m2   % =0 ?????????????????????
%*********************************************************************
    
DIFF_m3 = DIFF_m2./dz;           % Compute the differences in concentration due to advection [g/m^3/d]
TLCDFF_m3 = TLCDFF_m2/dz(nl_soil);           % Total leaching in [g/m^3/d]


    
    