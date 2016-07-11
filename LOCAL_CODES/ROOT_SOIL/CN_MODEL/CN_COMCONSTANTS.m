function [kd, kl, kh] = ...
        COMCONSTANTS(ADD, Cb, Cl, Ch, fSd, fTd, rh, rr, phi, SWITCHES)

%=========================================================================
% This code solves the linear system of equations (1) (8) and (11) in 
% Porportato (2003) under steady state conditions.
%
%------------------------- Input Variables -------------------------------
%         ADD           % C added to that layer [gC/m^3/d] 
%         Cb            % Total C concetration in the biomass pool for that layer % [gC/m^3]
%         Cl            % Total C concetration in the litter pool for that layer % [gC/m^3]
%         Ch            % Total C concetration in the humus pool for that layer % [gC/m^3]
%         fds           % Nondimensional factor that describes soil
%                         moisture effects on decomposition
%         rh            % Isohumic coefficient that represents the fraction
%                         of decomposition litter that undergoes humification
%         rr            % Portion of decomposing carbon that is lost by
%                        respiration
%         phi           % Nondimensional factor that accounts for a possible reduction of the decomposition rate when the litte   
%
%------------------------- Output Variables ------------------------------
%         kd            % Value that defines death of microbial biomass
%         kl            % Value that defines the rate of decomposition of
%                       litter
%         kh            % Value that defines the rate of decomposition of
%                       humus
%=========================================================================
% solution of Ak = B 

if SWITCHES.CN_type == 1 

        % Defines matrix A
        A(1,1) = Cb; 
        A(1,2) = -phi*fds*Cb*Cl;
        A(1,3) = 0;

        A(2,1) = 0;
        A(2,2) = (rh*phi*fds*Cb*Cl);
        A(2,3) = -(phi*fds*Cb*Ch);

        A(3,1) = -Cb;
        A(3,2) = (1-rh-rr)*(phi*fds*Cb*Cl);
        A(3,3) = (1-rr)*(phi*fds*Cb*Ch);

        %Defines matrix B
        B(1,1) =-ADD;
        B(2,1) = 0;
        B(3,1) = 0;

        %Compute the vector
        x=A^(-1)*B;

        kd=x(1);
        kl=x(2);
        kh=x(3);

elseif SWITCHES.CN_type == 0

    % Defines matrix A
    A(1,1) = -phi*fSd*fTd*Cb*Cl; 
    A(1,2) = Cb;

    A(2,1) = phi*fSd*fTd*Cb*Cl*(1-rr); 
    A(2,2) = -Cb;


    %Defines matrix B
    B(1,1) =-ADD;
    B(2,1) = 0;

    %Compute the vector
    x=A^(-1)*B;

    kl=x(1);
    kd=x(2);
    kh =nan;

end
