function [VARIABLES] = SOILHEAT (Hg, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS)   
                    
 
%=========================================================================
%  This function implements the numerical solution for soil heat transport
%  described in Oleson et al (2004) for the Community Land Model (CLM)
%  Top BC: Hg is heat flux into top of soil column
%  Bottom BC: zero heat flux at bottom of soil column
%
%  The following equation is solved:
%
%   dT    d        dT  
% C -- = -- ( TK * -- )
%   dt   dz        dz  
%   
%  Written By: Darren Drewry and Modified by Juan Quijano
%
%------------------------- Input Variables -------------------------------
%       Hg              % [W / m^2] Soil heat flux at top of soil column        
%       VARIABLES       % VARIABLES structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%       CONSTANTS       % CONSTANTS structure
%------------------------- Output Variables ------------------------------
%       VARIABLES       % VARIABLES structure
% 
%========================================================================
%                    
% LOCAL VARIABLES:
%   amx             "a" left off diagonal of tridiagonal matrix
%   bmx             "b" diagonal column for tridiagonal matrix
%   cmx             "c" right off diagonal tridiagonal matrix
%   rmx             "r" forcing term of tridiagonal matrix
%=========================================================================


%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************                    
% PARAMS
Tf = PARAMS.Soil.Tf - 273.15;                                              % [ºC] Frezing temperature of water 
nl_soil = PARAMS.nl_soil;                                                  % [] Number of layers in the soil 
alph = PARAMS.Soil.alphCN;                                                 % Cranck Nicholson for Heat Solution  
HC_liq = PARAMS.Soil.HC_liq;                                               % [J / kg / K] Heat Capacity Liquid Water  
HC_ice = PARAMS.Soil.HC_ice;                                               % [J / kg / K] Heat Capacity Ice 
rho_liq = PARAMS.Soil.rho_liq;                                             % [kg / m^3] Density Liquid Water 
rho_ice = PARAMS.Soil.rho_ice;                                             % [kg / m^3] Density Ice Water 
TK_liq = PARAMS.Soil.TK_liq;                                               % [W / m / K] Thermal Conductivity liquid water 
TK_ice = PARAMS.Soil.TK_ice;                                               % [W / m / K] Thermal Conductivity Ice water 

% CONSTANTS
grav =CONSTANTS.grav;                                                      % [m / s^2] Gravity Acceleration 
dt =CONSTANTS.dtime;                                                       % [s] Time Step 

% VERTSTRUC
psi0 = VERTSTRUC.psi0;                                                     % [mm] Minimum Soil Suction at Saturation 
bsw=VERTSTRUC.bsw;                                                         % [] Clapp-Hornberger "b" parameter 
HKsat = VERTSTRUC.HKsat;                                                   % [mm/s] Hydraulic Conductivity at Saturation 
znode = VERTSTRUC.zns;                                                     % [m] Node Soil Layer Depth  
zlayer = VERTSTRUC.zhs;                                                    % [m] Interface Soil Layer Depth
dz = VERTSTRUC.dzs;                                                        % [m] Soil Layer Thickness 
TK_dry = VERTSTRUC.TK_dry;                                                 % [W / m / K] Thermal Conductivity of Dry Soil 
TK_sol = VERTSTRUC.TK_sol;                                                 % [W / m / K] Thermal Conductivity of Mineral Solids 
HC_sol = VERTSTRUC.HC_sol;                                                 % [J / m^3 / K] Heat Capacity of Soils Solids 
poros = VERTSTRUC.porsl;                                                   % [] Porosity of Soils 

%VARIABLES
volliq = VARIABLES.SOIL.volliq;                                            % [] Volumetric Liquid Water Content  
volice = VARIABLES.SOIL.volice;                                            % [] Volumetric Ice Water Content 
Ts = VARIABLES.SOIL.Ts;                                                    % [C] Temperature Soils  

%Compute wice and wliq
wice = volice * rho_ice;                                                   % [kg / m^2] Ice Density per unit area 
wliq = volliq * rho_liq;                                                   % [kg / m^2] Liquid Water Density per unit area


% SOIL SURFACE ENERGY EXCHANGE
%*******
% FIX THIS WHEN CANOPY MODEL COUPLED TO SOIL !!!!
% HEAT FLUX FROM ATMOSPHERE 
%   (Hg and dH_dT
dHg_dT=0;
%*******

Ts_prev = Ts;

% VOLUMETRIC HEAT CAPACITY [J / m^3 / K]  (6.73) Vries (1963) 
%    cpv = HC_sol.*(1-poros) + ((volliq * rho_liq)./dz).*HC_liq + ((volice * rho_ice)./dz).*HC_ice;% - 10^6;
     cpv = HC_sol.*(1-poros) + ((volliq * rho_liq)).*HC_liq + ((volice * rho_ice)).*HC_ice - 0*10^6;
   
    
% THERMAL CONDUCTIVITIES AT NODES 
    nfrinds = find(Ts>=Tf); % Indices of non-frozen layers
    frinds = find(Ts<Tf);   % Indices of frozen layers
    
    % Saturated Thermal Conductivity [W / m / K]
        TK_sat(nfrinds) = (TK_sol(nfrinds).^(1-poros(nfrinds))) .* (TK_liq.^(poros(nfrinds)));
        TK_sat(frinds) = (TK_sol(frinds).^(1-poros(frinds))) .* (TK_liq.^(poros(frinds))) .* ...
                         (TK_ice.^(poros(frinds)-volliq(frinds)));
        TK_sat = TK_sat(:);
        
    % (6.64)
        Sr = (volliq + volice) ./ poros;
        if (max(Sr) > 1)
            disp('ERROR in SoilHeatTransport: Sr > 1');
            keyboard;
        end
        
    % Kersten Number (6.63)
        Ke(nfrinds) = log10(Sr(nfrinds))+1;
        Ke(frinds) = Sr(frinds);
        Ke(find(Ke<0)) = 0; 
        Ke = Ke(:);
        %ke = ones(12,1)*0.8;
    
    % Soil Thermal Conductivity [W / m / K] (6.58)
        TKsoil = Ke .* TK_sat + (1-Ke).*TK_dry;
        inds = find(Sr <= 10^-7);
        TKsoil(inds) = TK_dry(inds);
    
    
% THERMAL CONDUCTIVITIES AT LAYER INTERFACES (6.11)
    num = TKsoil(1:nl_soil-1) .* TKsoil(2:nl_soil) .* (znode(2:nl_soil) - znode(1:nl_soil-1));
    denom = TKsoil(1:nl_soil-1) .* (znode(2:nl_soil) - zlayer(1:nl_soil-1)) + ...
            TKsoil(2:nl_soil) .* (zlayer(1:nl_soil-1) - znode(1:nl_soil-1));
    TKsoil_h(1:nl_soil-1) = num ./ denom;
    TKsoil_h(nl_soil) = 0;
    TKsoil_h = TKsoil_h(:);
    

% SET UP TRI-DIAGONAL SOLUTION
    fact = dt ./ cpv ./ dz;
    zdiff = diff(znode);
    
    % TOP NODE:
        aa(1) = 0;
        bb(1) = 1 + fact(1) * ((1-alph)*TKsoil_h(1)/(znode(2)-znode(1)) - dHg_dT);
        cc(1) = -(1-alph) * fact(1) * TKsoil_h(1) / (znode(2)-znode(1));
        rr(1) = Ts(1) + fact(1) * (Hg - dHg_dT*Ts(1) - alph * TKsoil_h(1) * (Ts(1)-Ts(2))/(znode(2)-znode(1)));

    % INTERNAL NODES:
        inds = [2:nl_soil-1];
        t1 = TKsoil_h(inds-1) ./ (znode(inds) - znode(inds-1));
        t2 = TKsoil_h(inds) ./ (znode(inds+1) - znode(inds));
        aa(inds) = -(1-alph) * fact(inds) .* t1;
        bb(inds) = 1 + (1-alph) * fact(inds) .* (t1 + t2);
        cc(inds) = - (1-alph) * fact(inds) .* t2;
        
        Fi = -t2 .* (Ts(inds) - Ts(inds+1));
        Fim1 = -t1 .* (Ts(inds-1) - Ts(inds));
        rr(inds) = Ts(inds) + alph * fact(inds) .* (Fi - Fim1);

    % BOTTOM NODE:
        aa(nl_soil) = -(1-alph) * fact(nl_soil) * TKsoil_h(nl_soil-1) / (znode(nl_soil) - znode(nl_soil-1));
        bb(nl_soil) = 1 + (1-alph) * fact(nl_soil) * TKsoil_h(nl_soil-1) / (znode(nl_soil) - znode(nl_soil-1));
        cc(nl_soil) = 0;
        rr(nl_soil) = Ts(nl_soil) + alph * fact(nl_soil) * TKsoil_h(nl_soil-1) / ...
                    (znode(nl_soil) - znode(nl_soil-1)) * (Ts(nl_soil-1) - Ts(nl_soil));

                
% CALL TRI-DIAGONAL SOLVER
    Ts_new = TRIDIAG(nl_soil, aa, bb, cc, rr);
    Ts_new = Ts_new(:);
    
% COMPUTE THE HEAT FLUXES
    Gnew = -((Ts_new(2:nl_soil)-Ts_new(1:nl_soil-1)).*TKsoil_h(1:nl_soil-1))./(znode(2:nl_soil)-znode(1:nl_soil-1));

% CORRECT VALUES 

%[Ts_new] = correctheat (Ts_new, Ts_prev, Tf, TKsoil_h, cpv, volliq, volice,...
%            rho_liq, rho_ice, bsw, grav, psi0, znode, dz, dt, nl_soil, poros, alph,...
%            Hg, wice, wliq);
    
% STORE IN VARIABLES
   VARIABLES.SOIL.Gnew = [Gnew ; 0];                                        % [W / m^2] Heat Fluxes in the Soil
   VARIABLES.SOIL.Ts = Ts_new;                                              % [C] Soil Temperature
   VARIABLES.SOIL.cpv = cpv;                                                % [J / m^3 / K] Volumetric Heat Capacity of Soil
   VARIABLES.SOIL.TKsoil = TKsoil;                                          % [W / m / K] Thermal Conductivity Soil
   VARIABLES.SOIL.TKsol = TK_sol;                                           % [W / m / K] Thermal Conductivity of Mineral Solids
   