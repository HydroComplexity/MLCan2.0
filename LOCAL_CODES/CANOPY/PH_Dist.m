function [Vz] = PH_Dist( VERTSTRUC, PARAMS )

%=========================================================================
% Vertical distribution of photosynthetic capacity from (Leuning et al., 1995 --> Eqn. 9)
%
% Written By: Darren Drewry
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%------------------------- Output Variables ------------------------------
%       Vz              % Vertical distribution of photosynthetic capacity
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VERTSTRUC
    LAI = VERTSTRUC.LAIz;                                                  % [m^2 leaf/ m^2 ground] LAI for all layers 

    % PARMAS
    kn = PARAMS.Photosyn.kn_canopy;                                        % [] Paramater to compute the vertical Distribution of Photosynthetic Capacity
%*************************************************************************

    LAIcum = cumsum ( flipdim(LAI(:),1) );

    LAIcum_top = flipdim(LAIcum(:), 1);

    Vz = exp(-kn * LAIcum_top);
    
    Vz = [Vz(2:end); 1];