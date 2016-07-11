function [UP_nit_active, Fix_Nit] = CN_passive_fixation_nuptake_Nit (nl_soil, nspecies, Nit_DEM, SUP_nit, Nit_ku, Nit, NFr);

% for Nit
%     for i=1:nl_soil
%         Nit_CaseEqn31(i)=Nit_DEM(i)-SUP_nit(i);
%         if Nit_CaseEqn31(i) < 0 || Nit_CaseEqn31(i) == 0
%             UP_nit_active(i) = 0;  % [g/m2/d]
%             Fix_Nit(i) = 0;  % [g/m2/d]
%         end
%         if Nit_ku(i) < Nit_CaseEqn31(i)
%             rup=(Nit_DEM(i)-SUP_nit(i)-SUP_nit(i)*NFr)/(Nit_ku(i)*Nit(i)*(1+NFr));
%             if rup < 0
%                 rup=0;
%             end
%             if rup > 1
%                 rup=1;
%             end
%             if rup == nan
%                 rup=0;
%             end
%             
%             UP_nit_active(i) = (Nit_ku(i)*Nit(i))*rup;
%             Fix_Nit(i) = min((Nit_ku(i)*Nit(i)*rup+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-Nit_ku(i)*Nit(i)*rup);
%         end
%         if Nit_ku(i) > Nit_CaseEqn31(i) && Nit_CaseEqn31(i) >0
%             rup=(Nit_DEM(i)-SUP_nit(i)-SUP_nit(i)*NFr)/((Nit_DEM(i)-SUP_nit(i))*(1+NFr));
%             if rup < 0
%                 rup=0;
%             end
%             if rup > 1
%                 rup=1;
%             end
%             if rup == nan
%                 rup=0;
%             end
%             
%             UP_nit_active(i) = (Nit_DEM(i)-SUP_nit(i))*rup;
%             Fix_Nit(i) = min(((Nit_DEM(i)-SUP_nit(i))*rup+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-(Nit_DEM(i)-SUP_nit(i))*rup);
%         end
%     end 

for j=1:nspecies
    for i=1:nl_soil
        Nit_CaseEqn31(i,j)=Nit_DEM(i,j)-SUP_nit(i,j);
        if Nit_CaseEqn31(i,j) < 0 || Nit_CaseEqn31(i,j) == 0
            UP_nit_active(i,j) = 0;  % [g/m2/d]
            Fix_Nit(i,j) = 0;  % [g/m2/d]
        end
        if Nit_ku(i,j) < Nit_CaseEqn31(i,j)
            rup=(Nit_DEM(i,j)-SUP_nit(i,j)-SUP_nit(i,j)*NFr(j))/(Nit_ku(i,j)*Nit(i)*(1+NFr(j)));
            if rup < 0
                rup=0;
            end
            if rup > 1
                rup=1;
            end
            if rup == nan
                rup=0;
            end
            
            UP_nit_active(i,j) = (Nit_ku(i,j)*Nit(i))*rup;
            Fix_Nit(i,j) = min((Nit_ku(i,j)*Nit(i)*rup+SUP_nit(i,j))*NFr(j),Nit_DEM(i,j)-SUP_nit(i,j)-Nit_ku(i,j)*Nit(i)*rup);
        end
        if Nit_ku(i,j) > Nit_CaseEqn31(i,j) && Nit_CaseEqn31(i,j) >0
            rup=(Nit_DEM(i,j)-SUP_nit(i,j)-SUP_nit(i,j)*NFr(j))/((Nit_DEM(i,j)-SUP_nit(i,j))*(1+NFr(j)));
            if rup < 0
                rup=0;
            end
            if rup > 1
                rup=1;
            end
            if rup == nan
                rup=0;
            end
            
            UP_nit_active(i,j) = (Nit_DEM(i,j)-SUP_nit(i,j))*rup;
            Fix_Nit(i,j) = min(((Nit_DEM(i,j)-SUP_nit(i,j))*rup+SUP_nit(i,j))*NFr(j),Nit_DEM(i,j)-SUP_nit(i,j)-(Nit_DEM(i,j)-SUP_nit(i,j))*rup);
        end
    end
end