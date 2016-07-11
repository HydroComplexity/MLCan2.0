
% compute PHI

function [phi, PHI, MIN_net, IMM_net, MIN_gross, IMM_gross, Nreg, DECl] = CN_computephi (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, ADD, CNa) 

Cl = VARIABLES.Cl;%       Cl = carbon concentration in litter pool [gC / m^3]
Ch = VARIABLES.Ch;%       Ch = carbon concentration in humus pool [gC / m^3]
Cb = VARIABLES.Cb;%       Cb = carbon concentration in biomass pool [gC / m^3]
Amm = VARIABLES.Amm;%       Amm = ammonium concentration in soil [gN / m^3]
Nit = VARIABLES.Nit;%       Nit = nitrate concentration in soil [gN / m^3]
kl = PARAMS.CN.kl;%       kl = rate of decomposition of the litter pool [m^3 d / gC]
kh = PARAMS.CN.kh;%       kh = rate of decomposition of the humus pool [m^3 d / gC]
koae = PARAMS.CN.koae;% Organic Assimilation Efficiency parameter
rr = PARAMS.CN.rr;%       rr = fraction of decomposed organic carbon that goes to respiration [-]
CNl = VARIABLES.CNl;%       CNl = carbon/nitrogen ratio of litter pool [gC / gN]
CNb = PARAMS.CN.CNb;%       CNb = carbon/nitrogen ratio of biomass pool [gC / gN]
CNh = PARAMS.CN.CNh;%       CNh = carbon/nitrogen ratio of humus pool [gC / gN]
if SWITCHES.CN.Bioturbation;
    nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
else
    nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
end
ku_Amm = PARAMS.CN.ku_Amm;
ku_Nit = PARAMS.CN.ku_Nit;

% Decomposition rate of litter
ratel = phi .* fTd .* fSd .* kl .* Cb;
DECl = ratel .* Cl;        

for ii = 1:nl_soil   
        % Eqn. 18. Compute PHI
        phi(ii) = 1;
        if SWITCHES.CN_type ==1 
            PHI(ii) = phi*fTd(ii)*fSd(ii)*Cb(ii)*(kh(ii)*Ch(ii)*(1./CNh - (1-rr)./CNb(ii)) + kl(ii)*Cl(ii)*(1./CNl(ii) - rh(ii)./CNh - (1-rh(ii)-rr)./CNb(ii)));            
        elseif SWITCHES.CN_type == 2      
            PHI(ii) = phi(ii)*fTd(ii)*fSd(ii)*Cb(ii)*(kl(ii)*Cl(ii)*(1./CNl(ii) - (1-rr-rr*rh(ii))./CNb(ii)));
        else
            PHI(ii) = phi(ii)*fTd(ii)*fSd(ii)*Cb(ii)*kl(ii)*Cl(ii)*(koae./CNl(ii) - (1-rr)./CNb(ii));                
        end
        if (PHI(ii)>=0)   % NET MINERALIZATION
                        
            IMM_Amm(ii) = 0;
            IMM_Nit(ii) = 0;
            
            MIN_net(ii) = PHI(ii) + (1-koae)*DECl(ii)/CNl(ii);
            IMM_net(ii) = 0;            
            MIN_gross(ii) = PHI(ii) + (1-koae)*DECl(ii)/CNl(ii);
            IMM_gross(ii) = 0;
            
            Nreg(ii) = 1;
            
        else        % NET IMMOBILIZATION        
            
            MIN_net(ii) = 0;
            
            IMMmax = fSd(ii)*(ku_Amm*Amm(ii) + ku_Nit*Nit(ii));
            
            if (abs(PHI(ii))<IMMmax) % UNRESTRICTED IMMOBILIZATION
                
                phi(ii) = 1;                
                MIN_net(ii) = max(PHI(ii) + (1-koae)*DECl(ii)/CNl(ii),0);
                IMM_net(ii) = -min(PHI(ii) + (1-koae)*DECl(ii)/CNl(ii),0);            
                MIN_gross(ii) = (1-koae)*DECl(ii)/CNl(ii);
                IMM_gross(ii) = abs(PHI(ii));

                if MIN_net(ii) > 0 
                    Nreg(ii) = 2;
                else  
                    Nreg(ii) = 4;
                end
                
            else            % RESTRICTED IMMOBILIZATION
                
                num = IMMmax;
                if SWITCHES.CN_type == 1
                    denom = fTd(ii)*fSd(ii)*Cb(ii)*(kh(ii)*Ch(ii)*(1/CNh - (1-rr)/CNb(ii)) + kl(ii)*Cl(ii)*(1/CNl(ii) - rh(ii)/CNh - (1-rh(ii)-rr)/CNb(ii)));                    
                elseif SWITCHES.CN_type == 2
                    denom = fTd(ii)*fSd(ii)*Cb(ii)*kl(ii)*Cl(ii)*(1/CNl(ii) - (1-rr-rr*rh(ii))/CNb(ii));
                else
                    denom = fTd(ii)*fSd(ii)*Cb(ii)*kl(ii)*Cl(ii)*(koae/CNl(ii) - (1-rr)/CNb(ii));                    
                end
                phi(ii) = -num / denom;  
                
                PHI(ii) = -IMMmax;
                MIN_net(ii) = max(PHI(ii) + (1-koae)*DECl(ii)/CNl(ii),0);
                IMM_net(ii) = -min(PHI(ii) + (1-koae)*DECl(ii)/CNl(ii),0);            
                MIN_gross(ii) = (1-koae)*DECl(ii)/CNl(ii);
                IMM_gross(ii) = PHI(ii);                                
                
                if MIN_net(ii) > 0 
                    Nreg(ii) = 3;
                else  
                    Nreg(ii) = 5;
                end
            end
        end
end    
PHI=PHI(:);
MIN_net=MIN_net(:); 
IMM_net=IMM_net(:);
MIN_gross=MIN_gross(:);
IMM_gross=IMM_gross(:);
Nreg=Nreg(:);


