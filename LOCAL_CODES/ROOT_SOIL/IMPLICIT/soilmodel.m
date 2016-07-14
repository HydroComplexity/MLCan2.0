% Dongkook
%function [dwat,psicom,kf,kl,fluxt,fluxb,flux_s,flux_sr,flux_sr_all,mberrormm,type, hor_drainage, hor_drainage_lay] ...
function [dwat,psicom,kf,kl,fluxt,fluxb,flux_s,flux_sr,flux_sr_all,mberrormm,type, hor_drainage, hor_drainage_lay,flux_Ss] ...
= soilmodel(nl_soil,dtime,thetas,...
                            pentry,bpar,ks,zsoi,dzsoi,zisoi,...
                            thetai, plants, ...
                            psiroot,krad,pthr,nspecies, type, n, alpha, thetar,vanGen)

%=========================================================================
% This function solves the implicit scheme for a given time step. In this
% function an iteratively solution is implemented
%
% Written by Juan Camilo Quijano, UIUC, 2008
%
%------------------------- Input Variables -------------------------------
%       nl_soil         % [] Number of layes
%       dtime           % [s] Time step 
%       thetas          % [] Soil moisture at saturation 
%       pentry          % [mm] minimum soil suction, i.e., soil potential at saturation 
%       bpar            % Clapp-Hornbereger "b" parameter
%       ks              % [m/s] Hydraulic conductivity at saturation
%       zsoi            % [mm] Vector with depth at nodes of each layer
%       dzsoi           % [mm] Vector with Layer thicknesses 
%       zisoi           % [mm] Vector with depth of interfaces between the layers  
%       thetai          % [] Volumetric water content in previous time step 
%       plants          % [1 or 0] 1. Simulation with Plants 0. Simulation with no Plants        
%       psiroot         % [mm] Root Water Potential
%       krad            % Root Radial conductivities [1/s]
%       pthr            % Rainfall that reach the soil (Not intercepted). Flux [mm/s] 
%       nspecies        % [] Number of species
%       type            % 1.OK, 2. Supersaturation First Layer, 3. Super Saturation in between layers, thus hor_drainage >0 
%------------------------- Output Variables ------------------------------
%       dwat            % [] Change in volumetric water content
%       psicom          % [mm] Soil matrix solution
%       kf              % [mm/s] Hydraulic conductivity at interface between soil layers
%       fluxt           % [mm/s] Flux at the top
%       fluxb           % [mm/s] Flux at the bottom
%       flux_s          % [mm/s] Fluxes between soil layers
%       flux_sr         % [mm/s] Fluxes between soil and roots at all layers 
%       flux_sr_all     % [mm/s] Fluxes between soil and roots at all layers 
%       mberrormm       % [mm] Mass Balance error in Soil-Root solution.
%       type            % 1.OK, 2. Supersaturation First Layer, 3. super Saturation in between layers, thus hor_drainage >0 
%       hor_drainage    % [mm/s] Horizontal drainage if super saturation occurs
%       hor_drainage_lay% [mm/s] Horizontal drainage if super saturation for all layers
%-------------------------------------------------------------------------                        
%  Create vector to compute the flux boundary condition

vec1=[1 10^(-1) 10^(-2) 10^(-3) 10^(-4)];
vec2=(vec1(1:4)+vec1(2:5))/2;
indf = (1:5)*2-1; 
vecf(indf) = vec1;
indf1 = (1:5)*2-1; 
indf2 = (1:4)*2;
vecf(indf1) = vec1;
vecf(indf2) = vec2;
%-------------------------------------------------------------------------                        

% Dongkook: Start
%vanGen=1; % if it is 0, Brooks Corey
lambda  = n-1;                  % [-]
m       = lambda./n;             % [-]

%alpha   = 0.01*(1/10);         % parameter related to the inverse of the air entry suction [1/cm] to [1/mm]
%n       = 2.0;                  % Pore-size distributions [-]
%thetar = 0.08;                 % Residual water content [-]


if vanGen == 1
    [psii, ki, kli, thetai] = characteristic_vanGenuchten(thetai,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,2,alpha,n,lambda,m,thetar);
else
    [psii, ki, kli, thetai] = characteristic(thetai,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,2);
end
% Dongkook: End

for i=1:1:length(vecf)  

        criteria=10000; criteriamin=10000;  
        iteration=0;
        while criteria>0.05*10^(-1);
            iteration=iteration+1; 
            if (iteration>20);
                break;
            end;   % Initalize the variables
            if iteration==1
                psiant=psii;thetaant=thetai;kant=ki;
                % Dongkook, Start
                %[Cant] = computec(psiant,thetas,nl_soil,bpar,pentry); 
                if vanGen == 1
                    [Cant] = computecS_vanGenuchten(psiant,thetaant,thetas,nl_soil,bpar,pentry,alpha,n,lambda,m,thetar);                    
                else
                    [Cant] = computecS(psiant,thetaant,thetas,nl_soil,bpar,pentry);
                end
                % Dongkook - End
            else
                psiant=psicom; 
                
                if vanGen == 1
                    [thetaant, kant, klant, psiant] = characteristic_vanGenuchten(psiant,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,1,alpha,n,lambda,m,thetar);                    
                else                    
                    [thetaant, kant, klant, psiant] = characteristic(psiant,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,1);
                end
                % Dongkook, Start
                %[Cant] = computec(psiant,thetas,nl_soil,bpar,pentry);      
                if vanGen == 1
                    [Cant] = computecS_vanGenuchten(psiant,thetaant,thetas,nl_soil,bpar,pentry,alpha,n,lambda,m,thetar);                    
                else
                    [Cant] = computecS(psiant,thetaant,thetas,nl_soil,bpar,pentry);
                end
                % Dongkook - End
            end
           
            if vanGen == 1
                Se = ((thetaant - thetar)./(thetas - thetar));
                
                fb=ks(nl_soil)*Se(nl_soil)^(1/2).*(1 - (1 - Se(nl_soil)^(1/m(nl_soil)))^m(nl_soil))^2;
            else                
                fb=ks(nl_soil)*(thetaant(nl_soil)/thetas(nl_soil))^(2*bpar(nl_soil)+3);
            end
            
            %ks(nmc).*((vars(nmc))./thetas(nmc)).^(2.*bpar(nmc)+3);  
%            ft=min([pthr ((eff_porosity(1)-thetaant(1))*dzsoi(1))/dtime ks(1)]);
            ft=min([pthr*vecf(i) ks(1)]); 
            
            
            % Dongkook, Start
            % Based on Le et al., 20015 to modify the soil moisture equation
            %[Am,KK,GG,CC,KKr] = matrices(nl_soil,Cant,kant,dzsoi,dtime,ft,fb,zsoi,krad,nspecies);%Compute matrices
            [Am,KK,GG,CC,KKr,CC_Ss] = matrices(nl_soil,Cant,kant,dzsoi,dtime,ft,fb,zsoi,krad,nspecies,thetas,thetaant);%Compute matrices
            % Dongkook, End
             
%             psicom=(Am+CC+KKr)^(-1)*(-GG-KK+CC*psii+KKr*psiroot);
             KKrsol2 = zeros(nl_soil,1);
             KKrsol1= sum(KKr,3);
             for ii=1:1:nspecies
                 KKrsol2 = KKr(:,:,ii)*psiroot(:,ii) + KKrsol2;
             end
             
             % Dongkook - Start
             % Using "^(-1)" in MATLAB is not a good way to solve it. 
             % Let's use \ instaed.          
             %psicom=(Am+CC+KKrsol1)^(-1)*(GG+KK+CC*psiant-(1/dtime)*(thetaant-thetai)+KKrsol2);                  
             psicom=(Am+CC+KKrsol1+CC_Ss)\(GG+KK+CC*psiant+CC_Ss*psii-(1/dtime)*(thetaant-thetai)+KKrsol2);
             if psicom(1) > 0
                stop =10 ;
             end
             % Dongkook - End
             
             if vanGen == 1
                 [thetaf, kf, kl, psicom] = characteristic_vanGenuchten(psicom,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,1,alpha,n,lambda,m,thetar);                 
             else
                 [thetaf, kf, kl, psicom] = characteristic(psicom,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,1);
             end
%            criteria=max(abs((thetaf-thetaant))./abs(thetaant));
            criteria=max(abs((psicom-psiant))./abs(psiant));
            
            if criteria < criteriamin
                psicoms=psicom;
                criteriamin=criteria;
            end
        end
        if iteration < 20
            break
        end        
end
if i > 1
     type(2) = 1;
end
if (iteration > 20)
   psicom=psicoms;
end

if vanGen == 1
    [thetaf, kf, kl, psicom] = characteristic_vanGenuchten(psicom,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,1,alpha,n,lambda,m,thetar);    
else
    [thetaf, kf, kl, psicom] = characteristic(psicom,thetas,ks,nl_soil,zsoi,dzsoi,zisoi,bpar,pentry,1);
end
% Dongkook
%[flux_s,flux_sr,flux_sr_all,fluxt,fluxb,dwat,mberror,mberrormm]= ...
%compflux(psicom,psiroot,kf,krad,ft,fb,dzsoi,zsoi,nl_soil,type,thetai,thetaf,dtime, nspecies);

[flux_s,flux_sr,flux_sr_all,fluxt,fluxb,dwat,flux_Ss,mberror,mberrormm]= ...
compflux(psicom,psiroot,kf,krad,ft,fb,dzsoi,zsoi,nl_soil,type,thetai,thetaf,CC_Ss, psii, dtime, nspecies);
% Dongkook - End

% If there is supersaturation of layer in the sooil column other than the
% one in the top the excess water it released from the soil column as horizontal
% drainage

% Dongkook, Start
% if (max(psicom-pentry)>=0)
%     type(3)=1;
%     [flux_s,flux_sr,flux_sr_all,fluxt,fluxb,dwat,mberror,mberrormm, hor_drainage,hor_drainage_lay, psicom]= ...
%     correct_saturation(psicom,psiroot,kf,krad,ft,fb,dzsoi,zsoi,nl_soil,type,thetai,...
%     thetaf,dtime, nspecies, pentry, thetas);
% else
% Dongkook - End    
    hor_drainage = 0;   
    hor_drainage_lay = zeros(nl_soil,1);
% Dongkook, Start
%end
% Dongkook - End    

