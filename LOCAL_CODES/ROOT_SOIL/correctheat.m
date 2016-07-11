
function [Ts_new] = correctheat (Ts_new, Ts_prev, Tf, TKsoil_h, cpv, volliq, volice,...
            rho_liq, rho_ice, bsw, grav, psi0, znode, dz, dt, nl_soil, porsl, alph,...
            Hg, wice, wliq)

if sum(Ts_new < 0)>0;
    stop = 43;
end

% This function corrects the energy balance in the soil for those cases
% where melting or freezing occurs.

Lf = 3.337e5;            % latent Heat of fusion in [J/kg] 

% Compute the minimum wliq in ice

nfrinds = Ts_prev>=Tf; % Indices of non-frozen layers
frinds = Ts_prev<Tf;   % Indices of frozen layers

wliqmax(frinds) = dz(frinds).*porsl(frinds).*((10^3)*Lf*(Tf-Ts_prev(frinds))./...
    (grav.*Ts_prev(frinds).*psi0(frinds))).^(-(1./bsw(frinds)));
wliqmax(nfrinds) = 0;

% there is a mistake in CLM tecnical document. Below is the correction 
% using Nui and Yang formulation eq 3.
wliqmax = wliqmax .*rho_liq;
wliqmax = wliqmax(:);

% COMPUTE THE EXCESS OF DEFICIT OOF ENERGY NEEDED IN EACH LAYER 
% First layer
Fp_2 = -(Ts_prev(2)-Ts_prev(1))/(znode(2)-znode(1))*TKsoil_h(1);
Fn_2 = -(Ts_new(2)-Ts_new(1))/(znode(2)-znode(1))*TKsoil_h(1);
H(1) =  -(alph*(Fp_2)+(1-alph)*(Fn_2)-Hg) - cpv(1)*dz(1)*(Tf-Ts_prev(1))/(dt);

% From 2 to N-1 layer
Fp_2 = -(Ts_prev(3:nl_soil)-Ts_prev(2:nl_soil-1))./(znode(3:nl_soil)-znode(2:nl_soil-1)).*TKsoil_h(2:nl_soil-1);
Fp_1 = -(Ts_prev(2:nl_soil-1)-Ts_prev(1:nl_soil-2))./(znode(2:nl_soil-1)-znode(1:nl_soil-2)).*TKsoil_h(1:nl_soil-2);

Fn_2 = -(Ts_new(3:nl_soil)-Ts_new(2:nl_soil-1))./(znode(3:nl_soil)-znode(2:nl_soil-1)).*TKsoil_h(2:nl_soil-1);
Fn_1 = -(Ts_new(2:nl_soil-1)-Ts_new(1:nl_soil-2))./(znode(2:nl_soil-1)-znode(1:nl_soil-2)).*TKsoil_h(1:nl_soil-2);

H(2:nl_soil-1) = -(alph*(Fp_2-Fp_1) + (1-alph)*(Fn_2-Fn_1))...
    - cpv(2:nl_soil-1).*dz(2:nl_soil-1).*(Tf-Ts_prev(2:nl_soil-1))/(dt);
% At layer N
Fp_1 = -(Ts_prev(nl_soil)-Ts_prev(nl_soil-1))/(znode(nl_soil)-znode(nl_soil-1))*TKsoil_h(nl_soil-1);
Fn_1 = -(Ts_new(nl_soil)-Ts_new(nl_soil-1))/(znode(nl_soil)-znode(nl_soil-1))*TKsoil_h(nl_soil-1);
Fp_2 = 0;
Fn_2 = 0;

H(nl_soil) = -(alph*(Fp_2-Fp_1) + (1-alph)*(Fn_2-Fn_1))...
    - cpv(nl_soil)*dz(nl_soil)*(Tf-Ts_prev(nl_soil))/(dt);
H = H(:);

% Compute Hm
Hm = H*dt./Lf;                              % [kg]

% COMPUTE THE NEW VALUES OF WICE
wice_new = wice; % Initialize with previous values

% compute the Melting, Ti_(n+1) >Tf, Hm > 0, wice > 0 
ind = min([Ts_new > Tf  Hm > 0  wice > 0],[],2);
wice_new(ind) = wice(ind) - Hm(ind);

% compute the Melting, Ti_(n+1) < Tf, Hm < 0, wliq > wliqmax 
ind = min([Ts_new < Tf  Hm <= 0  wliq>wliqmax],[],2);
wice_new(ind) = min([wice(ind)+wliq(ind)-wliqmax(ind) wice(ind)-Hm(ind)],[],2);


% Readjust liquid water
deltawice = wice_new - wice;
wliq_new = wliq - deltawice;


if max(wice_new) > 0
    stop =5;
end

% Compute again the energy needed
Hstar = H + Lf*(deltawice)/dt; 

% Compute the new temperature 
Ts_new = Tf + dt*Hstar./cpv./dz;








