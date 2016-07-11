function [res LE H LW Tl]  = LEB_residual(ee,dd,bb,Tl,Ta,ea,Rabs,Me)

%=========================================================================
% This Function Computes the Resiude of the leaf energy balance.
%
% ¡¡¡¡¡¡¡¡¡
% - This formulation is used when there is no good closure using the
% exponential apporoximation
% - It is applied only for the wet fractions, in dry fraction the 
% closure problem does not arise
%
% Written By: Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       ee              % Factor in Nikolov et al (1995)  model
%       dd              % Factor in Nikolov et al (1995)  model
%       bb              % Factor in Nikolov et al (1995)  model
%       Tl              % [C] Leaf temperature profile
%       Ta              % [C] Air Temperature 
%       ea              % [kPa] Vapor Pressure Air,
%       Rabs            % [W / m^2] Absorbed Radiation in the leaf
%       Me              % [W / m^2] - energy stored in biochemical reactions
%------------------------- Output Variables ------------------------------
%       res             % [W / m^2] Residual in the Energy Balance 
%       LE              % [W / m^2] LE from leaf 
%       H               % [W / m^2] H from leaf 
%       LW              % [W / m^2] LW radiation emitted from leaf
%       Tl              % [C] Leaf temperature 
% 
%========================================================================



% Compute the saturated vapor pressure at Tl
es =  0.6107 * exp( (17.38*Tl) ./ (Tl + 239) );
%es = 610.7*exp((17.38*Tl)/(Tl+239))/1000;

% Cmopute the LE term
LE = ee*(es-ea);

% Compute the H term
H = dd*(Tl-Ta);

% Compute the LW emitted
LW = bb*(Tl+273.16)^4;

% Compute the residue
res = Rabs - LE - H - LW - Me;

