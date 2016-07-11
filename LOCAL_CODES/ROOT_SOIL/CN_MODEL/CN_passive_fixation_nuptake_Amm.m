function [UP_amm_active, Fix_Amm] = CN_passive_fixation_nuptake_Amm (nl_soil, nspecies, Amm_DEM, SUP_amm, Amm_ku, Amm, NFr);

% % for Amm
% for i=1:nl_soil    
%     % for Amm
%     Amm_CaseEqn31(i)=Amm_DEM(i)-SUP_amm(i);
%     if Amm_CaseEqn31(i) < 0 || Amm_CaseEqn31(i) == 0
%         UP_amm_active(i) = 0;  % [g/m2/d]
%         Fix_Amm(i) = 0;  % [g/m2/d]
%     end
%     if Amm_ku(i) < Amm_CaseEqn31(i)
%         rup=(Amm_DEM(i)-SUP_amm(i)-SUP_amm(i)*NFr)/(Amm_ku(i)*Amm(i)*(1+NFr));
%         if rup < 0
%             rup=0;
%         end
%         if rup > 1
%             rup=1;
%         end
%         if rup == nan
%             rup=0;
%         end
%         
%         UP_amm_active(i) = (Amm_ku(i)*Amm(i))*rup;
%         Fix_Amm(i) = min((Amm_ku(i)*Amm(i)*rup+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-Amm_ku(i)*Amm(i)*rup);
%     end
%     if Amm_ku(i) > Amm_CaseEqn31(i) && Amm_CaseEqn31(i) >0
%         rup=(Amm_DEM(i)-SUP_amm(i)-SUP_amm(i)*NFr)/((Amm_DEM(i)-SUP_amm(i))*(1+NFr));
%         if rup < 0
%             rup=0;
%         end
%         if rup > 1
%             rup=1;
%         end
%         if rup == nan
%             rup=0;
%         end
%         
%         UP_amm_active(i) = (Amm_DEM(i)-SUP_amm(i))*rup;
%         Fix_Amm(i) = min(((Amm_DEM(i)-SUP_amm(i))*rup+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-(Amm_DEM(i)-SUP_amm(i))*rup);
%     end
% end

% for Amm
for j=1:nspecies
    for i=1:nl_soil
        Amm_CaseEqn31(i,j)=Amm_DEM(i,j)-SUP_amm(i,j);
        if Amm_CaseEqn31(i,j) < 0 || Amm_CaseEqn31(i,j) == 0
            UP_amm_active(i,j) = 0;  % [g/m2/d]
            Fix_Amm(i,j) = 0;  % [g/m2/d]
        end
        if Amm_ku(i,j) < Amm_CaseEqn31(i,j)
            rup=(Amm_DEM(i,j)-SUP_amm(i,j)-SUP_amm(i,j)*NFr(j))/(Amm_ku(i,j)*Amm(i)*(1+NFr(j)));
            if rup < 0
                rup=0;
            end
            if rup > 1
                rup=1;
            end
            if rup == nan
                rup=0;
            end
            
            UP_amm_active(i,j) = (Amm_ku(i,j)*Amm(i))*rup;
            Fix_Amm(i,j) = min((Amm_ku(i,j)*Amm(i)*rup+SUP_amm(i,j))*NFr(j),Amm_DEM(i,j)-SUP_amm(i,j)-Amm_ku(i,j)*Amm(i)*rup);
        end
        if Amm_ku(i,j) > Amm_CaseEqn31(i,j) && Amm_CaseEqn31(i,j) >0
            rup=(Amm_DEM(i,j)-SUP_amm(i,j)-SUP_amm(i,j)*NFr(j))/((Amm_DEM(i,j)-SUP_amm(i,j))*(1+NFr(j)));
            if rup < 0
                rup=0;
            end
            if rup > 1
                rup=1;
            end
            if rup == nan
                rup=0;
            end
            
            UP_amm_active(i,j) = (Amm_DEM(i,j)-SUP_amm(i,j))*rup;
            Fix_Amm(i,j) = min(((Amm_DEM(i,j)-SUP_amm(i,j))*rup+SUP_amm(i,j))*NFr(j),Amm_DEM(i,j)-SUP_amm(i,j)-(Amm_DEM(i,j)-SUP_amm(i,j))*rup);
        end
    end
end