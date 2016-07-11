


if (SWITCHES.soilCN_on && SWITCHES.CN.Bioturbation)
    Cl_store = NaN(nl_soil+1,N);
    Ch_store = NaN(nl_soil+1,N);
    Cb_store = NaN(nl_soil+1,N);
    Nl_store = NaN(nl_soil+1,N);
    Amm_store = NaN(nl_soil+1,N);
    Nit_store = NaN(nl_soil+1,N);
    Nitdd_store = NaN(nl_soil+1,N);    
    Ammdd_store = NaN(nl_soil+1,N);
    CNl_store = NaN(nl_soil+1,N);
        
    MIN_net_store = NaN(nl_soil+1,N);
    MIN_gross_store = NaN(nl_soil+1,N);
    IMM_net_store= NaN(nl_soil+1,N);
    IMM_gross_store = NaN(nl_soil+1,N);
    Nreg_store = NaN(nl_soil+1,N);
    DECl_store = NaN(nl_soil+1,N);
    PHI_store = NaN(nl_soil+1,N);
    phi_store = NaN(nl_soil+1,N);
    fSd_store = NaN(nl_soil+1,N);
    fTd_store = NaN(nl_soil+1,N);
    mberrorN_store = NaN(1,N);
    mberrorC_store = NaN(1,N); 
            
    %new

    litterthickness_store = NaN(1,N);
    Cl_change_store = NaN(nl_soil+1,N);
    Cflux_store=NaN(1,N);
    Nflux_store=NaN(1,N);    
    Cl_bio_change_store=NaN(nl_soil+1,N);
    Cl_bio_in_store=NaN(nl_soil+1,N);
    Cl_bio_out_store=NaN(nl_soil+1,N);
    bioNerror_store=NaN(1,N);
    bioCerror_store=NaN(1,N);
    dCb_dt_store = NaN(nl_soil+1,N);
    dCl_dt_store = NaN(nl_soil+1,N);
    ADD_store = NaN(nl_soil+1,N);
    ADD_bio_store = NaN(nl_soil+1,N);
    ADD_ex_store = NaN(nl_soil+1,N);
    OUT_store = NaN(nl_soil+1,N);
    ADD_net_store = NaN(nl_soil+1,N);
        

    % NITRATE
    ADV_nit_m2_store = NaN(nl_soil+1,N);                      %[g/m^2/d]
    ADV_nit_m3_store = NaN(nl_soil+1,N);                      %[g/m^3/d]
    ADVr_nit_m2_store = NaN(nl_soil+1,N);                    %[g/m^2/d]
    ADVr_nit_m3_store = NaN(nl_soil+1,N);                    %[g/m^3/d]
    ADVrdd_nit_m3_store = NaN(nl_soil+1,N);                %[g/m^3/d]
    DIFF_nit_m3_store = NaN(nl_soil+1,N);                    %[g/m^3/d]
    DIFF_nit_m2_store = NaN(nl_soil+1,N);                    %[g/m^2/d]
    LCH_nit_m3_store = NaN(nl_soil+1,N);       %[g/m^3/d]
    LCH_nit_m2_store = NaN(nl_soil+1,N);       %[g/m^2/d]
    DIFFr_nit_m2_store = NaN(nl_soil+1,N);                      %[g/m^2/d]
    DIFFr_nit_m3_store = NaN(nl_soil+1,N);                      %[g/m^3/d]
    DIFFrdd_nit_m3_store = NaN(nl_soil+1,N);                  %[g/m^3/d]
    UP_nit_m3_store = NaN(nl_soil+1,N);                    %[g/m^3/d]
    UP_nit_m2_store = NaN(nl_soil+1,N);                    %[g/m^2/d]
    UP_nit_all_m2_store = NaN(nl_soil+1, nspecies, N);            %[g/m^2/d]
    errormbnit_store = NaN(1,N); 
    deltanit_store = NaN(nl_soil+1,N);         %[g/m^3/d]
    deltanitdd_store = NaN(nl_soil+1,N);      %[g/m^3/d]
       
    % AMMONIUM
    ADV_amm_m2_store = NaN(nl_soil+1,N);                      %[g/m^2/d]
    ADV_amm_m3_store = NaN(nl_soil+1,N);                      %[g/m^3/d]
    ADVr_amm_m2_store = NaN(nl_soil+1,N);                    %[g/m^2/d]
    ADVr_amm_m3_store = NaN(nl_soil+1,N);                    %[g/m^3/d]
    ADVrdd_amm_m3_store = NaN(nl_soil+1,N);                %[g/m^3/d]
    DIFF_amm_m3_store = NaN(nl_soil+1,N);                    %[g/m^3/d]
    DIFF_amm_m2_store = NaN(nl_soil+1,N);                    %[g/m^2/d]
    LCH_amm_m3_store = NaN(nl_soil+1,N);       %[g/m^3/d]
    LCH_amm_m2_store = NaN(nl_soil+1,N);       %[g/m^2/d]
    DIFFr_amm_m2_store = NaN(nl_soil+1,N);                      %[g/m^2/d]
    DIFFr_amm_m3_store = NaN(nl_soil+1,N);                      %[g/m^3/d]
    DIFFrdd_amm_m3_store= NaN(nl_soil+1,N);                  %[g/m^3/d]
    UP_amm_m3_store = NaN(nl_soil+1,N);                    %[g/m^3/d]
    UP_amm_m2_store = NaN(nl_soil+1,N);                    %[g/m^2/d]
    UP_amm_all_m2_store = NaN(nl_soil+1, nspecies, N);            %[g/m^2/d]
    errormbamm_store = NaN(nl_soil+1,N); 
    deltaamm_store = NaN(nl_soil+1,N);          %[g/m^3/d]
    deltaammdd_store = NaN(nl_soil+1,N);      %[g/m^3/d]
    
     
    % CARBON MASS BALANCE    
    Cinput_store = NaN(1,N);
    CBoutput_store = NaN(1,N);
    DECout_store = NaN(1,N);
    dCb_m2_store = NaN(1,N);
    dCl_m2_store = NaN(1,N);        
    
    % NITROGEN MASS BALANCE        
    UP_N_m2_store = NaN(1,N);
    LCH_N_m2_store = NaN(1,N); 
    Ninput_store = NaN(1,N);
    NBoutput_store = NaN(1,N);
    dNb_m2_store = NaN(1,N);
    dNl_m2_store = NaN(1,N);
    dAmm_m2_store = NaN(1,N);
    dNit_m2_store = NaN(1,N);
    dAmmdd_m2_store = NaN(1,N);
    dNitdd_m2_store = NaN(1,N);
    
    % ROOT ARQUITECTURE
    RMI_store = NaN(nl_soil+1,nspecies,N);
    RLI_store = NaN(nl_soil+1,nspecies,N);	
    RAI_store = NaN(nl_soil+1,nspecies,N);
    RVI_store = NaN(nl_soil+1,nspecies,N);
    RAID_store = NaN(nl_soil+1,nspecies,N);
    RVID_store = NaN(nl_soil+1,nspecies,N);
    
    
    % Denitrification
    Denitrif_store  = NaN(nl_soil+1,N); 
    Total_N20_store = NaN(nl_soil+1,N); 
    N2O_nit_store   = NaN(nl_soil+1,N); 
    
    
elseif SWITCHES.soilCN_on
    Cl_store = NaN(nl_soil,N);
    Ch_store = NaN(nl_soil,N);
    Cb_store = NaN(nl_soil,N);
    Nl_store = NaN(nl_soil,N);
    Amm_store = NaN(nl_soil,N);
    Nit_store = NaN(nl_soil,N);
    Nitdd_store = NaN(nl_soil,N);                      %[g/m^3/d]deltaammdd  = 
    Ammdd_store = NaN(nl_soil+1,N);
    CNl_store = NaN(nl_soil,N);
    dCb_dt_store = NaN(nl_soil+1,N);
    dCl_dt_store = NaN(nl_soil+1,N);
        
    MIN_net_store = NaN(nl_soil,N);
    MIN_gross_store = NaN(nl_soil,N);
    IMM_net_store= NaN(nl_soil,N);
    IMM_gross_store = NaN(nl_soil,N);
    Nreg_store = NaN(nl_soil,N);
    DECl_store = NaN(nl_soil,N);
    PHI_store = NaN(nl_soil,N);
    phi_store = NaN(nl_soil,N);
    fSd_store = NaN(nl_soil,N);
    fTd_store = NaN(nl_soil,N);
    mberrorN_store = NaN(nl_soil,N);
    mberrorC_store = NaN(nl_soil,N);
    
    %new
    Cflux_store = NaN(1,N);
    litterthickness_store = NaN(1,N);    
    
    ADD_store = NaN(nl_soil+1,N);
    
    
    % NITRATE
    ADV_nit_m2_store = NaN(nl_soil,N);                      %[g/m^2/d]
    ADV_nit_m3_store = NaN(nl_soil,N);                      %[g/m^3/d]
    ADVr_nit_m2_store = NaN(nl_soil,N);                    %[g/m^2/d]
    ADVr_nit_m3_store = NaN(nl_soil,N);                    %[g/m^3/d]
    ADVrdd_nit_m3_store = NaN(nl_soil,N);                %[g/m^3/d]
    DIFF_nit_m3_store = NaN(nl_soil,N);                    %[g/m^3/d]
    DIFF_nit_m2_store = NaN(nl_soil,N);                    %[g/m^2/d]
    LCH_nit_m3_store = NaN(nl_soil,N);       %[g/m^3/d]
    LCH_nit_m2_store = NaN(nl_soil,N);       %[g/m^2/d]
    DIFFr_nit_m2_store = NaN(nl_soil,N);                      %[g/m^2/d]
    DIFFr_nit_m3_store = NaN(nl_soil,N);                      %[g/m^3/d]
    DIFFrdd_nit_m3_store = NaN(nl_soil,N);                  %[g/m^3/d]
    UP_nit_m3_store = NaN(nl_soil,N);                    %[g/m^3/d]
    UP_nit_m2_store = NaN(nl_soil,N);                    %[g/m^2/d]
    UP_nit_all_m2_store = NaN(nl_soil, nspecies, N);            %[g/m^2/d]
    errormbnit_store = NaN(1,N); 
    deltanit_store = NaN(nl_soil,N);         %[g/m^3/d]
    deltanitdd_store = NaN(nl_soil,N);      %[g/m^3/d]


    % AMMONIUM
    ADV_amm_m2_store = NaN(nl_soil,N);                      %[g/m^2/d]
    ADV_amm_m3_store = NaN(nl_soil,N);                      %[g/m^3/d]
    ADVr_amm_m2_store = NaN(nl_soil,N);                    %[g/m^2/d]
    ADVr_amm_m3_store = NaN(nl_soil,N);                    %[g/m^3/d]
    ADVrdd_amm_m3_store = NaN(nl_soil+1,N);                %[g/m^3/d]
    DIFF_amm_m3_store = NaN(nl_soil,N);                    %[g/m^3/d]
    DIFF_amm_m2_store = NaN(nl_soil,N);                    %[g/m^2/d]
    LCH_amm_m3_store = NaN(nl_soil,N);       %[g/m^3/d]
    LCH_amm_m2_store = NaN(nl_soil,N);       %[g/m^2/d]
    DIFFr_amm_m2_store = NaN(nl_soil,N);                      %[g/m^2/d]
    DIFFr_amm_m3_store = NaN(nl_soil,N);                      %[g/m^3/d]
    DIFFrdd_amm_m3_store= NaN(nl_soil,N);                  %[g/m^3/d]
    UP_amm_m3_store = NaN(nl_soil,N);                    %[g/m^3/d]
    UP_amm_m2_store = NaN(nl_soil,N);                    %[g/m^2/d]
    UP_amm_all_m2_store = NaN(nl_soil, nspecies, N);            %[g/m^2/d]
    errormbamm_store = NaN(nl_soil,N); 
    deltaamm_store = NaN(nl_soil,N);          %[g/m^3/d]
    deltaammdd_store = NaN(nl_soil,N);      %[g/m^3/d]
    
        
    UP_N_m2_store = NaN(nl_soil,N);
    LCH_N_m2_store = NaN(nl_soil,N);    

    % ROOT ARQUITECTURE
    RMI_store = NaN(nl_soil,nspecies,N);
    RLI_store = NaN(nl_soil,nspecies,N);	
    RAI_store = NaN(nl_soil,nspecies,N);
    RVI_store = NaN(nl_soil,nspecies,N);
    RAID_store = NaN(nl_soil,nspecies,N);
    RVID_store = NaN(nl_soil,nspecies,N);
end