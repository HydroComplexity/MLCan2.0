%function [A,KK,GG,CC,KKr] = matrices(Ne,Ce,Ke,dz,dt,ft,fb,zsoi,krad,nspecies)
function [A,KK,GG,CC,KKr,CC_Ss] = matrices(Ne,Ce,Ke,dz,dt,ft,fb,zsoi,krad,nspecies,thetas,thetaant)

%=========================================================================
% This code computes the matrices that compose the linear system of the 
% fully implicit scheme. The matrices computed in this code are used to solve 
% richards equation with a flux boundary condition at the top and a flux
% boundary condition at the bottom. 
%
% Written by Juan Camilo Quijano, UIUC, 2008
%
%------------------------- Input Variables -------------------------------
%       Ne              % [] Number of layes
%       Ce              % Value of parameter C for each layer. Size (Ne x 1)
%       Ke              % Hydraulic conductivity in all the boundaries between layers
%                         Size (Ne-1 x 1)  [mm/s]
%       dz              % [mm] minimum soil suction, i.e., soil potential at saturation 
%       dt              % [s] Time step 
%       ft              % [mm/s] Flux Top boundary condition.
%       fb              % [mm/s] Flux boundary condition at the botoom. 
%       zsoi            % [mm] Vector with depth at nodes of each layer
%       krad            % Radial conductivities [1/s]
%       nspecies        % [] Number of species
%       thetas          % [] Soil moisture at saturation
%       thetaant        % [] Soil moisture
%------------------------- Output Variables ------------------------------
% This code gives as output the matrices A, KK, GG CC according to 
% equation   P'(A+CC)=(GG+KK+CC*P-(thn-th)/dt)
% where P' is matric pressure in new iteration, P is matric pressure in old iteration
% thn is soil moisture in previous iteration th is soil moisture in the
% previous time step and dt is the time of each time step.
%       A               % Matrix AA
%       KK              % Matrix KK
%       GG              % Matrix GG
%       CC              % Matrix CC
%       KKr             % Matrix KKr
%-------------------------------------------------------------------------                        

% Dongkook, Start
% Based on Le et al., 2015 to modify the soil moisture equation
Ss = 5e-4 / 1000 * ones(length(zsoi),1);  % Speficif storage coefficient [1/m] to [1/mm]. 
% Dongkook, End


%  Create vector to compute the flux boundary condition
% PREALLOCATE MEMORY
KK = nan(Ne,1);
KKr = nan(Ne,Ne,nspecies);

[den]=createden(Ne,zsoi);
%FIRST VALUES
%Matrix A
A(1,1)=(1/(dz(1)))*(Ke(1)/den(1));
A(1,2)=-(1/(dz(1)))*(Ke(1)/den(1));

A(Ne,Ne)=(1/(dz(Ne)))*(Ke(Ne-1)/den(Ne-1));
A(Ne,Ne-1)=-(1/(dz(Ne)))*(Ke(Ne-1)/den(Ne-1));

%Vector CC
CC(1,1)=(Ce(1)/dt); 
CC(Ne,Ne)=(Ce(Ne)/dt);

% Dongkook, Start
CC_Ss(1,1)=((Ss(1)/thetas(1)*thetaant(1))/dt); 
CC_Ss(Ne,Ne)=((Ss(Ne)/thetas(Ne)*thetaant(Ne))/dt);
% Dongkook, End

%Vector GG
GG=zeros(Ne,1);
GG(1)=ft/(dz(1));
GG(Ne)=-fb/(dz(Ne));

%Vector KK
KK(1)=-(1/(dz(1)))*(Ke(1)); 
KK(Ne,1)=(1/(dz(Ne)))*(Ke(Ne-1));

%THE REST FILLED WITH A LOOP
% The loop is row by row
for i=2:1:Ne-1;
    %Matrix A
    A(i,i-1)=-(1/(dz(i)))*(Ke(i-1)/den(i-1));     
    A(i,i+1)=-(1/(dz(i)))*(Ke(i)/den(i));   
    A(i,i)=(1/(dz(i)))*(Ke(i-1)/den(i-1)+Ke(i)/den(i));  
    %Vector KK
    KK(i,1)=-(1/(dz(i)))*(Ke(i)-Ke(i-1)); 
    
    %Vector CC
    CC(i,i)=(Ce(i))/(dt);
    
    % Dongkook, Start
    CC_Ss(i,i)=(Ss(i)/thetas(i)*thetaant(i))/(dt);
    % Dongkook, End
end

for i=1:1:nspecies
    kr=krad(:,i);
    kr=kr./dz';
    KKr(:,:,i) = diag(kr);
end
