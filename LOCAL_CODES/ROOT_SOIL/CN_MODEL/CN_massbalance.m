%*************************************************************************

if SWITCHES.CN.Bioturbation
    % CHECK MASS BALANCE OF C
    dz = [VARIABLES.SOIL.litterthickness ; VERTSTRUC.dzsmm/1000];%       dz = grid depth [m]
    Cinput = sum(ADD.*dz); %[grC/m2/d]  Input from Litter and Bioturbation
    CBoutput = sum(OUT.*dz);  %[grC/m2/d]  Output from Bioturbation
    DECout = sum(DECl.*dz*rr);   %[grC/m2/d] Respiration
    dCb_m2 = sum(dCb_dt.*dz);   %[grC/m2/d] change in Microbial Biomass 
    dCl_m2 = sum(dCl_dt.*dz);   %[grC/m2/d] change in organic matter
    mberrorC = Cinput - CBoutput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]
    
    % CHECK MASS BALANCE OF N
    Ninput = sum(ADD./CNa.*dz); %[gr/m2/d]
    NBoutput = sum(OUT./CNo.*dz); %[gr/m2/d]    
    dAmm_m2 = sum(dAmm_dt.*dz); %[gr/m2/d]
    dNit_m2 = sum(dNit_dt.*dz); %[gr/m2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        dAmmdd_m2 = sum(dAmmdd_dt.*RVID); %[gr/m2/d]
        dNitdd_m2 = sum(dNitdd_dt.*RVID); %[gr/m2/d]
    end
    dNb_m2 = sum(dNb_dt.*dz);   %[gr/m2/d]
    dNl_m2 = sum(dNl_dt.*dz);   %[gr/m2/d]
    inputN = Ninput - NBoutput - UP_N_m2 - LCH_N_m2;       %[gr/m2/d]
    
    if SWITCHES.CN.NupRootBiomass == 1
        outputN = dNb_m2 + dNl_m2 + dAmm_m2 + dNit_m2 + dAmmdd_m2 + dNitdd_m2; %[gr/m2/d]
    elseif SWITCHES.CN.NupRootBiomass == 0
        outputN = dNb_m2 + dNl_m2 + dAmm_m2 + dNit_m2; %[gr/m2/d]
    end
    mberrorN = inputN - outputN; %[gr/m2/d]
    % Dongkook - Edit
    if SWITCHES.CN_type == 1
        dCh_dt=dCh_dt;
        DECh=DECh;
        DEChout = sum(DECh.*dz.*rr);       %[grC/m2/d]
        dCh_m2 = sum(dCh_dt.*dz);   %[grC/m2/d]
        %mberrorC = Linput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]
        mberrorC = mberrorC - DEChout - dCh_m2;
        
        dNh_m2 = sum(dNh_dt.*dz);   %[gr/m2/d]
        % N Drainage
        %Drain_N_m2 = Drain_Nit./(200*200)+Drain_Amm./(200*200);
        mberrorN=mberrorN-dNh_m2;%-Drain_N_m2;
    end
    if SWITCHES.CN.Denitrification == 1
        Denitrif_m2 = sum(Denitrif.*dz); %[gr/m2/d]
        N2O_nit_m2 = sum(N2O_nit.*dz);
        mberrorN = mberrorN - Denitrif_m2 - N2O_nit_m2;
        
    end
    % Dongkook - Edit End
else
    % CHECK MASS BALANCE OF C
    Linput = sum(ADD.*(dz_mm/1000)); %[grC/m2/d]
    DECout = sum(DECl.*(dz_mm/1000)*rr);       %[grC/m2/d]
    dCb_m2 = sum(dCb_dt.*(dz_mm/1000));   %[grC/m2/d]
    dCl_m2 = sum(dCl_dt.*(dz_mm/1000));   %[grC/m2/d]
    mberrorC = Linput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]
    % CHECK MASS BALANCE OF N
    Linput = sum(ADD./CNa.*(dz_mm/1000)); %[gr/m2/d]
    dAmm_m2 = sum(dAmm_dt.*(dz_mm/1000)); %[gr/m2/d]
    dNit_m2 = sum(dNit_dt.*(dz_mm/1000)); %[gr/m2/d]
    dNb_m2 = sum(dNb_dt.*(dz_mm/1000));   %[gr/m2/d]
    dNl_m2 = sum(dNl_dt.*(dz_mm/1000));   %[gr/m2/d]
    inputN = Linput - UP_N_m2 - LCH_N_m2;       %[gr/m2/d]
    outputN = dAmm_m2 + dNit_m2 + dNb_m2 + dNl_m2; %[gr/m2/d]
    mberrorN = inputN - outputN; %[gr/m2/d]
end