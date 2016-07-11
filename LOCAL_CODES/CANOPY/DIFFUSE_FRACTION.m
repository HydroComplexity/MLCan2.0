function [fdiff] = DIFFUSE_FRACTION (zendeg, doy, SWdn)
%
%
%=========================================================================
%   Implements the algorithm of Spitters et al (1986) to calculate the
%   diffuse fraction of short wave radiation from incident global
%   shortwave, and the theoretical incident extra-terrestrial shortwave
% 
% Written By: Darren Drewry
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       zendeg           % [deg] zenith angle of sun
%       doy              % [] Day of year
%       SWdn             % [W / m^2 / s] Shortwave Radiation Down
%------------------------- Output Variables ------------------------------
%       fdiff            % [] Fraction of diffuse in incoming SW radiation  
%========================================================================              
%*************************************************************************

    dinds = find(zendeg<90);
    ninds = find(zendeg>=90);

    Scs = 1370; % [J/m^2/s] Solar Constant
    solelev = (90 - zendeg)*pi/180;
    So = Scs.*(1+0.033*cos((360*doy/365)*pi/180)) .* sin(solelev);

    So(ninds) = NaN;

    Sg_o_So = SWdn./So;

    % Fraction Diffuse
    sinbeta = sin(solelev);
    R = 0.847 - 1.61*sinbeta + 1.04*sinbeta.^2;
    K = (1.47 - R)./1.66;

    i1 = find(Sg_o_So <= 0.22);
    fdiff(i1) = 1;
    i2 = find(Sg_o_So > 0.22 & Sg_o_So <= 0.35);
    fdiff(i2) = 1 - 6.4*(Sg_o_So(i2) - 0.22).^2;
    i3 = find(Sg_o_So > 0.35 & Sg_o_So <= K);
    fdiff(i3) = 1.47 - 1.66*Sg_o_So(i3);
    i4 = find(Sg_o_So > K);
    fdiff(i4) = R(i4);

    