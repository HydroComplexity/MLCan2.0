% store data

% NUTRIENTS VARIABLES
    Cl_store(:,tt)=VARIABLES.Cl;       %[gC/m3]
    Ch_store(:,tt)=VARIABLES.Ch;       %[gC/m3]
    Cb_store(:,tt)=VARIABLES.Cb;       %[gC/m3]
    Nl_store(:,tt)=VARIABLES.Nl;        %[gN/m3] 
    Amm_store(:,tt)=VARIABLES.Amm;        %[gN/m3]
    Nit_store(:,tt)=VARIABLES.Nit;         %[gN/m3]
    Nitdd_store(:,tt) = VARIABLES.Nitdd;                      %[g/m^3/d]deltaammdd   
    Ammdd_store(:,tt) = VARIABLES.Ammdd;                      %[g/m^3/d]deltaammdd  =     
    CNl_store(:,tt)=VARIABLES.CNl;        %[gC/gN]    
    dCl_dt_store(:,tt) = VARIABLES.dCl_dt;  %[gC/m3/d]
    dCb_dt_store(:,tt) = VARIABLES.dCb_dt;  %[gC/m3/d]
    
    
    MIN_net_store(:,tt)=VARIABLES.MIN_net;     %[gN/m3/d]
    MIN_gross_store(:,tt)=VARIABLES.MIN_gross;     %[gN/m3/d]
    IMM_net_store(:,tt)=VARIABLES.IMM_net;      %[gN/m3/d]
    IMM_gross_store(:,tt)=VARIABLES.IMM_gross;     %[gN/m3/d]
    Nreg_store(:,tt)=VARIABLES.Nreg;
    DECl_store(:,tt)=VARIABLES.DECl;
    PHI_store(:,tt)=VARIABLES.PHI;
    phi_store(:,tt)=VARIABLES.phi; 
    fSd_store(:,tt)=VARIABLES.fSd;
    fTd_store(:,tt)=VARIABLES.fTd;
    mberrorN_store(:,tt)=VARIABLES.mberrorN;     %[gN/m3/d]
    mberrorC_store(:,tt)=VARIABLES.mberrorC;     %[gC/m3/d]    
    ADD_store(:,tt) = VARIABLES.ADD;            %[gr/m3/d]
    
    %new
    litterthickness_store(:,tt)=VARIABLES.SOIL.litterthickness;
    
    % NITRATE
    if SWITCHES.CN.NupRootBiomass == 1
        ADV_nit_m2_store(:,tt) = VARIABLES.ADV_nit_m2;                      %[g/m^2/d]
        ADV_nit_m3_store(:,tt) = VARIABLES.ADV_nit_m3;                      %[g/m^3/d]
        ADVr_nit_m2_store(:,tt) = VARIABLES.ADVr_nit_m2;                    %[g/m^2/d]
        ADVr_nit_m3_store(:,tt) = VARIABLES.ADVr_nit_m3;                    %[g/m^3/d]
        ADVrdd_nit_m3_store(:,tt) = VARIABLES.ADVrdd_nit_m3;                %[g/m^3/d]
        DIFF_nit_m3_store(:,tt) = VARIABLES.DIFF_nit_m3;                    %[g/m^3/d]
        DIFF_nit_m2_store(:,tt) = VARIABLES.DIFF_nit_m2;                    %[g/m^2/d]
    end
    LCH_nit_m3_store(:,tt) = VARIABLES.LCH_nit_m3;       %[g/m^3/d]
    LCH_nit_m2_store(:,tt) = VARIABLES.LCH_nit_m2;       %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        DIFFr_nit_m2_store(:,tt) = VARIABLES.DIFFr_nit_m2;                      %[g/m^2/d]
        DIFFr_nit_m3_store(:,tt) = VARIABLES.DIFFr_nit_m3;                      %[g/m^3/d]
        DIFFrdd_nit_m3_store(:,tt) = VARIABLES.DIFFrdd_nit_m3;                  %[g/m^3/d]
    end
    UP_nit_m3_store(:,tt) = VARIABLES.UP_nit_m3;                    %[g/m^3/d]
    UP_nit_m2_store(:,tt) = VARIABLES.UP_nit_m2;                    %[g/m^2/d]
    UP_nit_all_m2_store(:,:,tt) = VARIABLES.UP_nit_all_m2;            %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        errormbnit_store(:,tt) = VARIABLES.errormbnit;
        deltanit_store(:,tt) = VARIABLES.deltanit;          %[g/m^3/d]
        deltanitdd_store(:,tt) = VARIABLES.deltanitdd;      %[g/m^3/d]
    end
    
    % Dongkook Woo - Edit
    STORAGE.UP_nit_m2_store=UP_nit_all_m2_store;
    % Dongkook Woo - Edit End

    % AMMONIUM
    if SWITCHES.CN.NupRootBiomass == 1
        ADV_amm_m2_store(:,tt) = VARIABLES.ADV_amm_m2;                      %[g/m^2/d]
        ADV_amm_m3_store(:,tt) = VARIABLES.ADV_amm_m3;                      %[g/m^3/d]
        ADVr_amm_m2_store(:,tt) = VARIABLES.ADVr_amm_m2;                    %[g/m^2/d]
        ADVr_amm_m3_store(:,tt) = VARIABLES.ADVr_amm_m3;                    %[g/m^3/d]
        ADVrdd_amm_m3_store(:,tt) = VARIABLES.ADVrdd_amm_m3;                %[g/m^3/d]
        DIFF_amm_m3_store(:,tt) = VARIABLES.DIFF_amm_m3;                    %[g/m^3/d]
        DIFF_amm_m2_store(:,tt) = VARIABLES.DIFF_amm_m2;                    %[g/m^2/d]
    end
    LCH_amm_m3_store(:,tt) = VARIABLES.LCH_amm_m3;       %[g/m^3/d]
    LCH_amm_m2_store(:,tt) = VARIABLES.LCH_amm_m2;       %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        DIFFr_amm_m2_store(:,tt) = VARIABLES.DIFFr_amm_m2;                      %[g/m^2/d]
        DIFFr_amm_m3_store(:,tt) = VARIABLES.DIFFr_amm_m3;                      %[g/m^3/d]
        DIFFrdd_amm_m3_store(:,tt) = VARIABLES.DIFFrdd_amm_m3;                  %[g/m^3/d]
    end
    UP_amm_m3_store(:,tt) = VARIABLES.UP_amm_m3;                    %[g/m^3/d]
    UP_amm_m2_store(:,tt) = VARIABLES.UP_amm_m2;                    %[g/m^2/d]
    UP_amm_all_m2_store(:,:,tt) = VARIABLES.UP_amm_all_m2;            %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        errormbamm_store(:,tt) = VARIABLES.errormbamm;
        deltaamm_store(:,tt) = VARIABLES.deltaamm;          %[g/m^3/d]
        deltaammdd_store(:,tt) = VARIABLES.deltaammdd;      %[g/m^3/d]
    end
    
    % Dongkook Woo - Edit
    STORAGE.UP_amm_m2_store=UP_amm_all_m2_store;
    % Dongkook Woo - Edit End
    
    % CARBON MASS BALANCE    
    Cinput_store(tt) = VARIABLES.Cinput;
    CBoutput_store(tt) = VARIABLES.CBoutput;
    DECout_store(tt) = VARIABLES.DECout;
    dCb_m2_store(tt) = VARIABLES.dCb_m2;
    dCl_m2_store(tt) = VARIABLES.dCl_m2;    
        
    % NITROGEN MASS BALANCE   
    UP_N_m2_store(tt) = VARIABLES.UP_N_m2;
    LCH_N_m2_store(tt) = VARIABLES.LCH_N_m2;  
    Ninput_store(tt) = VARIABLES.Ninput;
    NBoutput_store(tt) = VARIABLES.NBoutput;
    dNb_m2_store(tt) = VARIABLES.dNb_m2;
    dNl_m2_store(tt) = VARIABLES.dNl_m2;
    dAmm_m2_store(tt) = VARIABLES.dAmm_m2;
    dNit_m2_store(tt) = VARIABLES.dNit_m2;
    if SWITCHES.CN.NupRootBiomass == 1
        dAmmdd_m2_store(tt) = VARIABLES.dAmmdd_m2;
        dNitdd_m2_store(tt) = VARIABLES.dNitdd_m2;
    end
        
    % ROOT ARQUITECTURE
    if SWITCHES.CN.NupRootBiomass == 1
        RMI_store(:,:,tt) = VARIABLES.RMI;
        RLI_store(:,:,tt) = VARIABLES.RLI;
        RAI_store(:,:,tt) = VARIABLES.RAI;
        RVI_store(:,:,tt) = VARIABLES.RVI;
        RAID_store(:,:,tt) = VARIABLES.RAID;
        RVID_store(:,:,tt) = VARIABLES.RVID;
    end
    
    % Dongkook Woo - Edit
    % Denitrification
    Denitrif_store(:,tt)  = VARIABLES.Denitrif;
    Total_N20_store(:,tt) = VARIABLES.Total_N20;
    N2O_nit_store(:,tt)   = VARIABLES.N2O_nit;
    
    Nitrif_store(:,tt) = VARIABLES.Nitrif;
    DECh_store(:,tt) = VARIABLES.DECh;
    % Dongkook Woo - End
%     
if SWITCHES.CN.Bioturbation   
    Cflux_store(:,tt) = VARIABLES.Cbiorate;
    Nflux_store(:,tt) = VARIABLES.Nbiorate;    
    Cl_bio_change_store(:,tt) = VARIABLES.Cl_bio_change;
    Cl_bio_in_store(:,tt) = VARIABLES.Cl_bio_in;   % [gr/m3/dtime]
    Cl_bio_out_store(:,tt) = VARIABLES.Cl_bio_out; % [gr/m3/dtime]
    bioCerror_store(:,tt) = VARIABLES.bioCerror;
    bioNerror_store(:,tt) = VARIABLES.bioNerror;
    ADD_bio_store(:,tt) = VARIABLES.ADD_bio;     %[gr/m3/d]
    ADD_ex_store(:,tt) = VARIABLES.ADD_ex;       %[gr/m3/d]
    OUT_store(:,tt) = VARIABLES.OUT;       % Only from bioturbation[gr/m3/d]
    ADD_net_store(:,tt) = VARIABLES.ADD_net;       %[gr/m3/d]
end
%     
