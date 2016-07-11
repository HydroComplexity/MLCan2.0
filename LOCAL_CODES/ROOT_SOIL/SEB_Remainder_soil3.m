function [remain, LEs, Gsl, Gs, Tlb] = SEB_Remainder_soil3(Tl, Ts1,...
                                   zsl, TK_litter, dzs1, TC1, Rabs)

                              

%=========================================================================
% Solve surface energy balance below a snow pack. (No Latent Heat)
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       Tl              % [C] Temperature in the Snow-Litter pack 
%       Tsl             % [C] Temperature in the Soil fist layer pack 
%       zsl             % [m] Depth snow pack 
%       TK_litter       % [W / m / K] Thermal Conductivity Snow Pack 
%       dzs1            % [m] Thickness first soil layer 
%       TC1             % [W / m / K] Thermal Conductivity Soil 
%       Rabs            % [W/m2] Energy Absorption Soil
%------------------------- Output Variables ------------------------------
%       remain          % [W/m2] Energy Balance Error   
%       LE              % [W/m2] Latent heat from soil (should be zero)   
%       Gs1             % [W/m2] Ground heat flux into the interface   
%       Gs              % [W/m2] Ground Heat Flux into the Soil     
%       Tlb             % [C] Temperature at the interface  
% 
%========================================================================                               
                                                              
        
    Tlb = (Ts1*TC1*zsl + Tl*TK_litter*dzs1 + Rabs/2*dzs1*zsl)/ ...
         (TC1*zsl + TK_litter*dzs1);
    
    % COMPUTE THE GROUND HEAT FLUX FROM THE LITTER TO THE SOIL LITTER
    % SURFACE
    LEs = 0;
    
    Gsl = (TK_litter)*(Tl - Tlb)/((zsl)/2); 

    Gs = (TC1)*(Tlb - Ts1)/((dzs1)/2); 
    
    remain = Gsl - LEs - Gs + Rabs;
    