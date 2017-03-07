      module eqn_state_m
c
        use e3_param_m
        use number_def_m
        use matdat_def_m
        use mmatpar_m
        use dgifinp_m
c
        implicit none
c
        real*8 :: rho_ref, p_ref, T_ref, alpha_P, beta_T, cv_liq
        real*8 :: rho_ref_s, p_ref_s, T_ref_s, alpha_P_s, cv_s
        real*8 :: bulkMod_s, shearMod_s 
c
      contains
c
      subroutine getthm6_ideal_gas
c
        mw    = mat_prop(mater,iprop_ideal_gas_mw, 1)
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = Ru/mw*1.0d3
        gamma1 = gamma - one
c
        rho = pres / (Rgas*T)
        ei  = T * Rgas / gamma1
c
      end subroutine getthm6_ideal_gas
c
      subroutine getthm7_ideal_gas
c
        call getthm6_ideal_gas
c
        h   = T * Rgas / gamma1 * gamma
        cp  = Rgas*gamma / gamma1
        alphaP = one / T
        betaT  = one / pres
        if (associated(cv)) cv  = Rgas / gamma1
        if (associated(gamb)) gamb = gamma1
        if (associated(c)) c =  sqrt( gamma * Rgas * T )
c
      end subroutine getthm7_ideal_gas
c
      subroutine getthm6_ideal_gas_mixture
c
        integer :: iel
        real*8, dimension(npro) :: mw,rgas,y
c
        y = vap_frac
        y = max(zero,y)
        y = min(one, y)
c
        mw  = y*MW_liquid + (one-y)*mat_prop(mater,iprop_ideal_gas_mw, 1)
c
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = Ru/mw*1.0d3
        gamma1 = gamma - one
c
        rho = pres / (Rgas*T)
        ei  = T * Rgas / gamma1
c
      end subroutine getthm6_ideal_gas_mixture
c
      subroutine getthm7_ideal_gas_mixture
c
        real*8, dimension(npro) :: mw,rgas,y
c
        y = vap_frac
        y = max(zero,y)
        y = min(one, y)
c
        mw  = y*MW_liquid + (one-y)*mat_prop(mater,iprop_ideal_gas_mw, 1)
c
        gamma = mat_prop(mater,iprop_ideal_gas_gamma,1)
        Rgas  = Ru/mw*1.0d3
        gamma1 = gamma - one
c
        rho = pres / (Rgas*T)
        ei  = T * Rgas / gamma1
c
        h   = T * Rgas / gamma1 * gamma
        cp  = Rgas*gamma / gamma1
        alphaP = one / T
        betaT  = one / pres
        if (associated(cv)) cv  = Rgas / gamma1
        if (associated(gamb)) gamb = gamma1
        if (associated(c)) c =  sqrt( gamma * Rgas * T )
c
      end subroutine getthm7_ideal_gas_mixture
c
      function rho_ideal_gas(p,R,T) result(rho)
        implicit none
        real*8 p,R,T,rho
        rho = p / (R*T)
      end function rho_ideal_gas
c
      subroutine getthm6_liquid_1
c
        rho_ref = mat_prop(mater,iprop_liquid_1_rho_ref,1)
        p_ref   = mat_prop(mater,iprop_liquid_1_p_ref,  1)
        T_ref   = mat_prop(mater,iprop_liquid_1_T_ref,  1)
        cv_liq  = mat_prop(mater,iprop_liquid_1_cv,     1)
        alpha_P = mat_prop(mater,iprop_liquid_1_alphaP, 1)
        beta_T  = mat_prop(mater,iprop_liquid_1_betaT,  1)
c
        rho = rho_ref * (one - alpha_P*(T-T_ref) + beta_T*(pres-P_ref))
        ei  = cv_liq*T
c
      end subroutine getthm6_liquid_1
c
      subroutine getthm7_liquid_1
c
        call getthm6_liquid_1
c
        h   = ei + pres/rho
        cp  = cv_liq
        alphaP = alpha_P
        betaT  = beta_T
        if (associated(cv)) cv  = cv_liq
        if (associated(gamb)) gamb = zero
c        c =  sqrt(one/(rho_ref*betaT))
        if (associated(c)) c =  sqrt(one/(rho*betaT))
c
      end subroutine getthm7_liquid_1
c
      subroutine getthm6_solid_1
c
        use e3_solid_m
        use number_def_m
        use propar_m 
        implicit none
c
        real*8, dimension(npro) :: temp_s
c
        rho_ref_s = mat_prop(mater,iprop_solid_1_rho_ref,1)
        p_ref_s   = mat_prop(mater,iprop_solid_1_p_ref,  1)
        T_ref_s   = mat_prop(mater,iprop_solid_1_T_ref,  1)
        cv_s     = mat_prop(mater,iprop_solid_1_cv,     1)
        alpha_P_s = mat_prop(mater,iprop_solid_1_alphaP, 1)
        bulkMod_s = mat_prop(mater,iprop_solid_1_bulkMod, 1)
        shearMod_s = mat_prop(mater,iprop_solid_1_shearMod, 1)
c
c        alphaP = alpha_P_s
c        betaT  = one /(bulkMod_s * Ja_def) ! double check here
c        rho = rho_ref_s * (one - alphaP*(T - T_ref_s) 
c     &                    + betaT*(pres - p_ref_s))
         stress_T_Mod = -three*bulkMod_s*alpha_P_s
         temp_s = stress_T_Mod *(T - T_ref_s) - bulkMod_s
     &          + (pres - p_ref_s) 
         rho = - (rho_ref_s * bulkMod_s)/(temp_s)
         alphaP = -(one/rho) *(rho_ref_s * bulkMod_s*stress_T_Mod)
     &          /(stress_T_Mod *(T - T_ref_s) - bulkMod_s)**two
         betaT = (one/rho) * (rho_ref_s * bulkMod_s)
     &         /( (pres - p_ref_s)-bulkMod_s )**two
c        ei  = ( cv_s - pres * alpha_P_s/rho)* T 
c    &        + (betaT * pres - alphaP * T)/rho * pres
        ei  = cv_s * T
c
      end subroutine getthm6_solid_1
c
      subroutine getthm7_solid_1
c
        use e3_solid_m
        implicit none
c
        call getthm6_solid_1
c
        h   = ei + pres/rho
        cv  = cv_s
        cp  = cv_s
        bulkMod = bulkMod_s
        shearMod = shearMod_s
c        c =  sqrt(one/(rho_ref*betaT))
C         c =  sqrt(one/(rho*betaT))
C         gamb = zero
c
c
      end subroutine getthm7_solid_1
c
      subroutine getthmif0
        use e3if_param_m
        call e3if_setparam_0
        call getthmif0_ptr
      end subroutine getthmif0
c
      subroutine getthmif1
        use e3if_param_m
        call e3if_setparam_1
        call getthmif1_ptr
      end subroutine getthmif1
c
      subroutine getthmif_solid_0
c
        use e3if_param_m
        use e3if_solid_data_m
        implicit none
c
        call e3if_setparam_solid_0
        call getthmif0_ptr            
c
      end subroutine getthmif_solid_0
c
      subroutine getthmif_solid_1
c
        use e3if_param_m
        use e3if_solid_data_m
        implicit none
c
        call e3if_setparam_solid_1
        call getthmif1_ptr  
c
      end subroutine getthmif_solid_1
c
      end module eqn_state_m
