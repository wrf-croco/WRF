!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_hs

CONTAINS

! GLmod ------------------------------------------------------------

! Ajout de la SST et humidite comme parametres de la routine
!    et en sortie GLW et GSW
!       GSW est relie a l'insolation...
   SUBROUTINE HSRAD(RTHRATEN,GLW,GSW,p8w,p_phy,pi_phy,dz8w,t_phy,          &
                     t8w, rho_phy, R_d,G,CP,dt,xlat,degrad,        &
                     QV, tsk,                           &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     )

!   SUBROUTINE HSRAD(RTHRATEN,p8w,p_phy,pi_phy,dz8w,t_phy,          &
!                     t8w, rho_phy, R_d,G,CP,dt,xlat,degrad,        &
!                     ids,ide, jds,jde, kds,kde,                    &
!                     ims,ime, jms,jme, kms,kme,                    &
!                     its,ite, jts,jte, kts,kte                     )





   IMPLICIT NONE

   INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte  

   REAL, INTENT(IN    )      ::        DEGRAD

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                              RTHRATEN
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(INOUT)  ::                                   GSW, &
                                                             GLW

   REAL, INTENT(IN   )   ::                   R_d,CP,G,dt

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                     dz8w, &
                                                             p8w, &
                                                           p_phy, &
                                                          pi_phy, &
                                                           t_phy, &
                                                             t8w, &
                                                         rho_phy  
   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN ) ::                                     xlat

   ! Cas Sec ou Humide
   LOGICAL :: drycase

   ! Parametres generaux
   INTEGER :: i,j,K,kte_,K_
   REAL :: dely,rcp,g0, r_terre  

   ! Parametres du Transfert Radiatif
   REAL :: sigma, nk

   ! Definition de l'epaisseur optique
   LOGICAL :: is_tau_front, is_tau_invert
   INTEGER :: create_sst
   REAL :: tau0pole,tau0eq,radflin,radfnlin
   REAL :: ytau, ltau, ysst, lsst, ssteq, sstpole
   LOGICAL :: use_tsk_tau, use_tsk_rad_bc
   REAL :: deltath_z

   ! Intermediaires pour le calcul
   REAL :: tau0, arg_phi, sst, tsfc
   REAL :: s8w
   REAL :: Dsw,Dp,Q_tend
   REAL :: n_tau,dtau
   REAL :: tau_strato,tau_tropo

   ! Flux montant et descendant
   ! Emissions des couches
   ! Epaisseurs optiques
   REAL, DIMENSION( kms:kme) :: F_up, F_down, B, tau

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT)  ::                   QV 

!   REAL, DIMENSION( ims:ime, jms:jme ),                           &
!         INTENT(INOUT)  ::                                    Q2 

   REAL, DIMENSION( ims:ime, jms:jme ),                           &
         INTENT(IN)  ::                                       tsk 



! A.Foussard
! We want to replace the relaxation with a simple
! LW radiative scheme, only forced by a SST field.
! This scheme give an equilibrium temperature different from
! HS'one.
! For details :
! A Gray-Radiation Aquaplanet Moist GCM. Part I:
! Static Stability and Eddy Scale
! 
! v5 :
! Pour avoir un champ atmospherique proche de celui des 
! simulations Bxx, uniquement a partir d'un champ de SST en
! equilibre radiatif avec l atmosphere
! Champ de sst plus frontal que pour Bxx
! Changement de flin pour ne pas avoir de destabilisation
! radiative de l'atmosphere.
! Egalement une epaisseur optique frontale (d'apres les calculs 
! de champ a l'equilibre, peu d'influence sur U et theta
! Reste une difference sur l'effet beta :
! Dans HS il n'est pas localise au contraire de Bxx
!
! v6 :
! On specifie une forme differente pour l epaisseur optique
! qui est obtenue par inversion du champ de temperature a 
! l'equilibre radiatif de Held et Suarez
!
! v7 : 
! Passer la SST comme argument de la fonction HSRAD
! Flux calcules sur les full level
!
! grey01 :
! Dans la troposphere : epaisseur optique
! qui est obtenue par inversion du champ de temperature a l 
! equilibre radiatif de Held et Suarez
! Dans la stratosphere : parametrisation de Frierson


!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level    F_down = 0
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level    T, dtau, dSW, RHTHRATEN, B
!         kms+2  ----- full level    tau, F_up, F_down, p8w, s8w
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level    
!         kms    ----- full level    F_up = sigma tsk**4
!
! The horizontal velocities used in the physics are unstaggered 
! relative to temperature/moisture variables. All predicted 
! variables are carried at half levels except w, which is at full 
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!==================================================================


   CALL nl_get_asiv_ssteq(1,ssteq)
   CALL nl_get_asiv_tau0eq(1,tau0eq)
   CALL nl_get_asiv_sstpole(1,sstpole)
   CALL nl_get_asiv_tau0pole(1,tau0pole)
   CALL nl_get_asiv_is_tau_front(1, is_tau_front)
   CALL nl_get_asiv_create_sst(1, create_sst)
   CALL nl_get_asiv_is_tau_invert(1, is_tau_invert)
   CALL nl_get_asiv_use_tsk_tau(1, use_tsk_tau)
   CALL nl_get_asiv_use_tsk_rad_bc(1, use_tsk_rad_bc)
   CALL nl_get_asiv_ltau(1,ltau)
   CALL nl_get_asiv_ytau(1,ytau)
   CALL nl_get_asiv_lsst(1,lsst)
   CALL nl_get_asiv_ysst(1,ysst)
   CALL nl_get_asiv_drycase(1,drycase)
   CALL nl_get_asiv_radflin(1,radflin)
   CALL nl_get_dy(1,dely)
   CALL nl_get_asiv_deltath_z(1,deltath_z)
 

   rcp      = R_d/CP
   g0       = 9.8
   n_tau    = 4.               ! exposant de la partie non lineaire
   sigma    = 5.6734E-8        ! Constante de Stefan W/m-2/K-4
   kte_     = MIN(kte,kde-1)
   nk       = 4. * rcp         ! exposant qui apparait
   radfnlin = 1. - radflin


   j_loop: DO J=jts,MIN(jte,jde-1)
   i_loop: DO I=its,MIN(ite,ide-1)

      ! CAS SEC
      LOOP_Q: DO K=kts,kte_
         ! Cas Sec : 
         ! mettre le rapport de melange a 0 sur les niveaux masse
         IF ( drycase ) THEN
           QV(I,K,J) = 0.
           ! RTHCUTEN(I,K,J) = 0.
           ! Q2(I,J) = 0.
         END IF
      ENDDO LOOP_Q
      

      ! Initialisation de l epaisseur optique
      IF ( is_tau_front ) THEN
         arg_phi = (j*dely/1000.0-ytau)/ltau
         tau0 = tau0eq - 0.5 * (tau0eq-tau0pole) * (1+tanh(arg_phi))
      ELSE
         ! epaisseur optique telle que dans Fr05
         arg_phi = (J*dely/1000.0-ytau)/4.0e3
         if(arg_phi>1.)arg_phi=1.
         if(arg_phi<-1.)arg_phi=-1.
         tau0=(tau0eq+tau0pole)/2.-(tau0eq-tau0pole)* &
                    sin(3.14159265*arg_phi/2.)/2.
      END IF


      IF ( use_tsk_tau ) THEN

         tsfc = (((1.0+tau0)/(2.0+tau0))**0.25) * tsk(I,J)

      ELSE            


      ! au cas ou...
         sst=293.

      ! Cas ou on specifie la SST afin de calculer l epaisseur optique
         arg_phi = (dely*float(J)/1000.0-ysst)/lsst

         IF ( create_sst.eq.0) THEN
            sst = 0.5*(ssteq+sstpole)

         ELSEIF(create_sst.eq.1)THEN
            arg_phi = (J*dely/1000.0-ysst)/4.0e3
            if(arg_phi>1.)arg_phi=1.
            if(arg_phi<-1.)arg_phi=-1.
            sst=(ssteq+sstpole)/2.-(ssteq-sstpole)* &
                    sin(3.14159265*arg_phi/2.)/2.

        ELSEIF(create_sst.eq.2)THEN
           sst= ssteq - 0.5*(ssteq - sstpole)*(1+tanh(arg_phi))

        END IF

         tsfc = (((1.0+tau0)/(2.0+tau0))**0.25) * sst

      END IF


! ------------------------------------------------------------
!                   CALCUL DES FLUX RADIATIFS


!   calcul de l'epaisseur optique et de l'emissivite de la couche

      LOOP_TAU: DO K=kts,kte_

         ! Niveau sigma : calcul sur les full levels (comme w)
         s8w=p8w(i,K,j)/p8w(i,1,j)        
      
         IF ( is_tau_invert ) THEN

            ! Epaisseur optique obtenue "par inversion" pour obtenir
            ! a l equilibre radiatif la meme structure verticale
            ! du champ de temperature que dans HS94 ou dans B05.
            ! Ici la stratification verticale est la meme a l equateur 
            ! ou au pole - pas le cas dans HS94
            tau_tropo = (1.0+tau0) * (s8w**nk) * ((1.0 - deltath_z*alog(s8w)/tsfc)**4) - 1.0
            tau_strato = tau0 * radflin * s8w

            ! Pas de temperature limite pour definir la tropopause,
            ! la condition tau > tau0*flin*sigma suffit a changer la stratification
            tau(K) = max(tau_strato, tau_tropo)

         ELSE

            ! Epaisseur optique telle que dans Fr05 avec la possibilite
            ! de faire varier le parametre radflin qui controle
            ! stabilite et stratification verticale
            tau(K) = tau0*(radflin*s8w + radfnlin*(s8w**n_tau))

         END IF 

         ! t_phy : Temperature absolue pour l'emissivite
         B(K) = sigma*(t_phy(i,k,j)**4.0)

      ENDDO LOOP_TAU



! ------------------------------------------------------------
!                 CONDITIONS AUX LIMITES

      ! Flux a la surface : 

      IF ( use_tsk_rad_bc ) THEN
! on utilise la sst reelle =>  effet radiatif des tourbillons
         F_up(kts) = sigma*(tsk(I,J)**4.0)
      ELSE
! sinon, on utilise une pseudo-sst
         F_up(kts) = sigma*(sst**4.0)
      END IF

      ! Flux entrant au sommet de l atmosphere
      F_down(kte_) = 0.    


! ------------------------------------------------------------
!               FLUX RADIATIFS HAUT ET BAS

      LOOP_FLUXES: DO K=kts,kte_-1

         ! Rayonnement vers le haut de l atmosphere 
         ! Les flux sont calcules sur la coordonnee verticale
         ! aux full levels W - grille staggered
         ! Besoin du flux montant aux half levels P : schema implicite
         dtau         = tau(K+1)-tau(K)    ! < 0 
         F_up(K+1)    = F_up(K) * (1.0 + 0.5*dtau) - B(K) * dtau
         F_up(K+1)    = F_up(K+1) / (1.0 - 0.5*dtau)

         ! Rayonnement vers la surface      
         ! Rayonnement integre de haut en bas : indices inverses
         K_           = kte_ + kts - K
         dtau         = tau(K_-1)-tau(K_) ! > 0 
         F_down(K_-1) = F_down(K_) * (1.0 - 0.5*dtau) + B(K_-1) * dtau
         F_down(K_-1) = F_down(K_-1) / (1.0 + 0.5*dtau)

      ENDDO LOOP_FLUXES


!       calcul du flux longwave arrivant a la surface
      GLW(I,J)=F_down(kts)
      GSW(I,J)=0.
! GLend


! ------------------------------------------------------------
!                     TENDANCE EN TEMPERATURE

      ! Bilan energetique sur chaque couche d'atmosphere
      ! Au niveau des half levels avec les flux calcules
      ! au full level : schema tres simple
      LOOP_RTHRATEN: DO K=kts,kte_-1


         Dsw       = F_up(K+1)-F_up(K)+F_down(K)-F_down(K+1)
         Dp        = p8w(i,K+1,j)-p8w(i,K,j)
         Q_tend        = g0/CP* Dsw/Dp

         ! On ajoute la valeur obtenue a la tendance calculee
         ! a partir de routine
         ! La tendance Qtend est en temperature ABSOLUE
         ! Terme en pi_phy : present dans ra_hs initial : 
         ! fonction d Exner, qui permet de passer de la 
         ! temperature absolue a la temperature potentielle
         ! Analyse dimensionelle : les wrf_rst 
         ! donnent RTHRATEN en K.s-1.Pa car la valeur 
         ! obtenue est multipliee par mu la pression totale de la
         ! colonne dans le coeur dynamique

         RTHRATEN(I,K,J)=RTHRATEN(I,K,J) + &
                         Q_tend / pi_phy(i,K,j)     

      ENDDO LOOP_RTHRATEN

   ENDDO i_loop
   ENDDO j_loop                                          





   END SUBROUTINE HSRAD


! GLend  ------------------------------------------------------------------

!====================================================================
   SUBROUTINE hsinit(RTHRATEN,restart,                              &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte                   )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           :: restart
   INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::        &
                                                          RTHRATEN


  INTEGER :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN
     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHRATEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

   ENDIF

   END SUBROUTINE hsinit

!====================================================================

END MODULE module_ra_hs



! GLmod --------------------------------------------------
!    on deplace tout ici...



!   INTEGER :: i,j,K,NK
!   real :: delty,delthez,p0,sec_p_d,sigb,kka,kks,kkf,rcp
!   real :: ttmp,teq,sig,sigterm,kkt,t_tend
!
! Newtonian relaxation scheme from Held and Suarez, Bull. Amer. Met.
! Soc., Vol. 75, No. 10., p1825-1830, 1994.  (box on page 1826)
! CEN and MIR  31-JUL-04
!
!   delty   = 60.0
!   delthez = 10.0
!   p0      = 100000.0
!   sec_p_d = 86400.
!   sigb    = 0.7
!   kka     = 1.0/40.0   ! units of per day
!   kks     = 0.25
!   kkf     = 1.0
!   rcp     = R_d/CP
!
!   j_loop: DO J=jts,MIN(jte,jde-1)
!   k_loop: DO K=kts,MIN(kte,kde-1)
!   i_loop: DO I=its,MIN(ite,ide-1)
!
!      ttmp = 315.0 - delty*(sin(xlat(i,j)*degrad))**2.0- &
!               delthez*alog(p_phy(i,k,j)/p0)*(cos(xlat(i,j)*degrad))**2.0
!               
!      teq=max(200.0,ttmp*(p_phy(i,k,j)/p0)**rcp)
!
!      sig=p_phy(i,k,j)/p8w(i,1,j)
!      sigterm=max(0.0,(sig-sigb)/(1.0-sigb))
!
!      kkt=kka+(kks-kka)*sigterm*(cos(xlat(i,j)*degrad))**4.0
!
!      t_tend=-kkt*(t_phy(i,k,j)-teq)/sec_p_d  ! t_tend in kinetic K/s
!
!      RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+t_tend/pi_phy(i,k,j)
!
!   ENDDO i_loop
!   ENDDO k_loop
!   ENDDO j_loop                                          
