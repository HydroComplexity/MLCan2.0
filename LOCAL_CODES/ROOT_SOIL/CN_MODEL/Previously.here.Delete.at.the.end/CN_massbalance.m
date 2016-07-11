%*************************************************************************

if SWITCHES.CN.Bioturbation
    % CHECK MASS BALANCE OF C
    dz = dz_mm/1000;
    Linput = sum(ADD.*dz); %[grC/m2/d]  Input from Litter and Bioturbation
    Boutput = sum(OUT.*dz);  %[grC/m2/d]  Output from Bioturbation
    DECout = sum(DECl.*dz*rr);   %[grC/m2/d] Respiration
    dCb_m2 = sum(dCb_dt.*dz);   %[grC/m2/d] change in Microbial Biomass 
    dCl_m2 = sum(dCl_dt.*dz);   %[grC/m2/d] change in organic matter
    mberrorC = Linput - Boutput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]
    
    % CHECK MASS BALANCE OF N
    Linput = sum(ADD./CNa.*dz); %[gr/m2/d]
    Boutput = sum(OUT./CNo.*dz); %[gr/m2/d]    
    dAmm_m2 = sum(dAmm_dt.*dz); %[gr/m2/d]
    dNit_m2 = sum(dNit_dt.*dz); %[gr/m2/d]
    dNb_m2 = sum(dNb_dt.*dz);   %[gr/m2/d]
    dNl_m2 = sum(dNl_dt.*dz);   %[gr/m2/d]
    inputN = Linput - Boutput - UP_N_m2 - LCH_N_m2;       %[gr/m2/d]
    outputN = dAmm_m2 + dNit_m2 + dNb_m2 + dNl_m2; %[gr/m2/d]
    mberrorN = inputN - outputN; %[gr/m2/d]
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