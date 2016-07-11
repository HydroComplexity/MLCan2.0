function [VERTSTRUC] ...
            = PEDOSOIL_PROPERTIES(PARAMS, VERTSTRUC, VARIABLES)
        
        
% This function computes Soil properties as a function of Organic 
% Matter and texture using pedotransfer functions described in 
% Saxton and Rawls 2006

% De reference structures
nl_soil = VERTSTRUC.nl_soil;

load pedotra_info;   % Load Soil Carbon and Bulk Densities in [g/m3]
clay = VERTSTRUC.clay.*ones(nl_soil,1);
sand = VERTSTRUC.sand.*ones(nl_soil,1);
Cl_soil = Cl(2:13);
rhod = rho_bulk';
OM_conc = Cl_soil*1.72;
OM = (OM_conc./rhod)*100;

% Check maximum OM equal to 8
indom = OM>8;
OM(indom) = 8;

P=PEDO_TRANSFER_FUNCTIONS(sand/100,clay/100,OM); % From Saxon and Rawls 2006 
                                                % sand and clay is in decimal
                                                % while OM is in %w look
                                                % figure 2 figure 3. It has
                                                % to be true to get numbers
                                                % in Table 3 
%******************************  TEST VALUES ******************************                                                
% Set Ksat in the organic horizon to 300 [mm/s] 
HKsat = P.Ksat/3600;   % from [mm/hr] to [mm/s]
HKsat(1) = 200/3600;                                               
% Test if psi0 is higher than zero seto to 0.0001 MPa
indpsi = P.hb > -0.0001;
P.hb(indpsi) = - 0.0001;                                                 
%**************************************************************************                                                
                                                
% ASSIGN
VERTSTRUC.porsl = P.Ss;    % POROSITY
VERTSTRUC.psi0 = P.hb/9.8066e-06;               % MINIMUM SOIL SUCTION = SOIL POTENTIAL AT SATURATION [mm]
VERTSTRUC.bsw = P.lam;                          % Clapp-Hornberger "b" parameter [-]
%VERTSTRUC.rhod = rhod;                          % BULK DENSITY [gr / m3]
VERTSTRUC.HKsat = HKsat;                  % HYDRAULIC CONDUCTIVITY AT SATURATION [mm / s] 
% VERTSTRUC.theta_dry = theta_dry;        
VERTSTRUC.eff_poros = VERTSTRUC.porsl;
                                                

                                                
                                           
       