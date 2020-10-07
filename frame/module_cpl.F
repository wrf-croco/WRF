MODULE module_cpl

   USE module_domain          , ONLY : domain, get_ijk_from_grid
   USE module_configure       , ONLY : grid_config_rec_type
   USE module_model_constants , ONLY : stbolt
   USE module_driver_constants, ONLY : max_domains, max_cplfld, max_extdomains
   USE module_cpl_oasis3 

   IMPLICIT NONE
   PRIVATE

   PUBLIC cpl_init
   PUBLIC cpl_set_dm_communicator
   PUBLIC cpl_defdomain
   PUBLIC cpl_settime
   PUBLIC cpl_snd
   PUBLIC cpl_rcv_sfcdrv
   PUBLIC cpl_store_input
   PUBLIC cpl_finalize
   PUBLIC cpl_abort

#ifdef key_cpp_oasis3
   LOGICAL     , PARAMETER, PUBLIC :: coupler_on = .TRUE.
   CHARACTER(5), PARAMETER         :: coupler_name = 'oasis'
#else
   LOGICAL     , PARAMETER, PUBLIC :: coupler_on = .FALSE.
   CHARACTER(4), PARAMETER         :: coupler_name = 'none'
#endif
   INTEGER :: nsecrun             ! current time in seconds since simulation restart
   INTEGER, PARAMETER :: charlen = 64
   CHARACTER(charlen), DIMENSION(max_domains,max_extdomains,max_cplfld) :: rcvname, sndname   ! coupling fields names for each nest

   CHARACTER(256) :: cltxt        ! messages or debug string
   INTEGER :: nlevdbg  = 1        ! verbosity level
   INTEGER :: nlevdbg2 = 10       ! verbosity level

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
   INCLUDE 'mpif.h'               ! only for MPI_COMM_NULL
#else
   INTEGER :: MPI_COMM_NULL = -1  ! define a fake (and not used) MPI_COMM_NULL, so it is compiling
#endif

CONTAINS

   SUBROUTINE cpl_init( kl_comm ) 
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_init  ***
      !!
      !! ** Purpose :   initialise coupling field names and WRF-coupler MPI communications
      !!--------------------------------------------------------------------
      INTEGER, INTENT(OUT) :: kl_comm       ! local MPI communicator of the model
      !
      INTEGER       :: jwrf,jext,jfld       ! local loop indicees
      CHARACTER( 3) :: clwrfdom, clextdom   ! d<domain>
      CHARACTER(16) :: clprefix             ! 'WRF_d??_EXT_d??_'
      !!--------------------------------------------------------------------

      ! coupling field name default definition
      rcvname(:,:,:) = 'not defined'
      sndname(:,:,:) = 'not defined'
      
      ! we could imagine to define rcvname and sndname through the namelist...
      ! define all possible coupling names with _d<domain> of WRF and the external model(s)
      DO jext = 1, max_extdomains
         
         WRITE(clextdom, fmt="('d',i2.2)") jext
         
         DO jwrf = 1, max_domains
            
            WRITE(clwrfdom, fmt="('d',i2.2)") jwrf          
            ! do not change following syntaxe as it is used in routines bellow 
            clprefix = 'WRF_'//clwrfdom//'_EXT_'//clextdom//'_' 
            
            !=========================================================================================================
            ! Variables that can be received by WRF :
            !=========================================================================================================
            
            ! surface driver: see subroutine cpl_rcv_sfcdrv (below) for details
            !--------------------------------------------------------------------------------------------------
            rcvname(jwrf,jext, 1) = clprefix//'SST'                  ! receive sea surface temperature (K)
            rcvname(jwrf,jext, 2) = clprefix//'UOCE'                 ! receive ocean zonal      surface current (m/s)
            rcvname(jwrf,jext, 3) = clprefix//'VOCE'                 ! receive ocean meridional surface current (m/s)
            rcvname(jwrf,jext, 4) = clprefix//'EOCE'                 ! receive ocean eastward   surface current (m/s)
            rcvname(jwrf,jext, 5) = clprefix//'NOCE'                 ! receive ocean northward  surface current (m/s)
            rcvname(jwrf,jext, 6) = clprefix//'CHA_COEF'             ! receive charnock coefficient
            
            !=========================================================================================================
            ! Variables that can be sent by WRF : see subroutine cpl_snd (below) for calculation details
            !=========================================================================================================
            
            ! mass fluxes (kg/m2/s = mm/s, positive downward)
            !--------------------------------------------------------------------------------------------------
            sndname(jwrf,jext, 1) = clprefix//'LIQUID_PRECIP'        ! send total liquid precipitation
            sndname(jwrf,jext, 2) = clprefix//'SOLID_PRECIP'         ! send total solid  precipitation
            sndname(jwrf,jext, 3) = clprefix//'TOTAL_EVAP'           ! send total evaporation: evap + sublimation
            sndname(jwrf,jext, 4) = clprefix//'EVAP-PRECIP'          ! send net fresh water budget: evap - total precipitation 
            
            ! heat fluxes (W.m-2, positive downward)
            !--------------------------------------------------------------------------------------------------
            sndname(jwrf,jext, 5) = clprefix//'SURF_NET_SOLAR'       ! send net surface shortwave heat flux
            sndname(jwrf,jext, 6) = clprefix//'SURF_NET_LONGWAVE'    ! send net surface  longwave heat flux
            sndname(jwrf,jext, 7) = clprefix//'SURF_LATENT'          ! send     surface    latent heat flux
            sndname(jwrf,jext, 8) = clprefix//'SURF_SENSIBLE'        ! send     surface  sensible heat flux
            sndname(jwrf,jext, 9) = clprefix//'SURF_NET_NON-SOLAR'   ! send net surface non-solar heat flux

            ! momentum fluxes (N m-2)
            !--------------------------------------------------------------------------------------------------
            sndname(jwrf,jext,10) = clprefix//'TAUX'                 ! send zonal      wind tress at atm-ocean interface
            sndname(jwrf,jext,11) = clprefix//'TAUY'                 ! send meridional wind tress at atm-ocean interface
            sndname(jwrf,jext,12) = clprefix//'TAUE'                 ! send Eastward   wind tress at atm-ocean interface
            sndname(jwrf,jext,13) = clprefix//'TAUN'                 ! send Northward  wind tress at atm-ocean interface
            sndname(jwrf,jext,14) = clprefix//'TAUMOD'               ! send module of  wind tress at atm-ocean interface

            ! wind speed for wave coupling (m s-1)
            !--------------------------------------------------------------------------------------------------
            sndname(jwrf,jext,15) = clprefix//'WND_U_01'             ! send the zonal wind speed at atmosphere-ocean interface
            sndname(jwrf,jext,16) = clprefix//'WND_V_01'             ! send the meridional wind speed at atmosphere-ocean interface
            sndname(jwrf,jext,17) = clprefix//'WND_E_01'             ! send the eastward wind speed at atmosphere-ocean interface
            sndname(jwrf,jext,18) = clprefix//'WND_N_01'             ! send the northward wind speed at atmosphere-ocean interface
            
         END DO
      END DO
      
      IF ( coupler_name == 'oasis' ) CALL cpl_oasis_init( kl_comm ) 
      
   END SUBROUTINE cpl_init
   

   SUBROUTINE cpl_set_dm_communicator( kdm_comm )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_initquilt  ***
      !!
      !! ** Purpose : provide the computing nodes communicator to the coupler
      !!--------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kdm_comm       ! MPI communicator between the computing nodes
      !!--------------------------------------------------------------------

      IF ( coupler_name == 'oasis' ) THEN 
         IF ( kdm_comm == MPI_COMM_NULL ) THEN
            CALL cpl_oasis_define( sndname, rcvname )   ! define io_quilting to OASIS
         ELSE
            CALL cpl_oasis_def_dmcomm( kdm_comm )       ! send the computing nodes communicator to OASIS
         END IF
      END IF

   END SUBROUTINE cpl_set_dm_communicator


   SUBROUTINE cpl_defdomain( grid )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_defdomain  ***
      !!
      !! ** Purpose : define each variable involved in the coupling and the grid partitioning
      !!--------------------------------------------------------------------
      TYPE(domain), INTENT(IN), POINTER ::   grid
      !
      INTEGER :: jwrf,jext,jfld          ! local loop indicees
      REAL    :: zmin,zmax               ! min/max of grid*cplmask
      INTEGER :: ips,ipe,jps,jpe,kps,kpe ! domain dimension on 1 processor
      INTEGER :: ims,ime,jms,jme,kms,kme ! memory domain dimension on 1 processor 
      INTEGER :: ids,ide,jds,jde,kds,kde ! domain dimension
      !!--------------------------------------------------------------------
#if (EM_CORE == 1)

      CALL get_ijk_from_grid( grid, ids, ide, jds, jde, kds, kde, &
         &                          ims, ime, jms, jme, kms, kme, &
         &                          ips, ipe, jps, jpe, kps, kpe  )

      ! first do some checks and prints. note that this could not be done in cpl_init
      ! which is called too early in the code

      ! some control prints on potential sent/received fields...
      CALL wrf_debug(nlevdbg, 'cpl_init: defined variables to be potentially received' )
      DO jfld = 1, max_cplfld
         DO jext = 1, grid%num_ext_model_couple_dom
            DO jwrf = 1, grid%max_dom
               IF( TRIM(sndname(jwrf,jext,jfld)) /= 'not defined' ) THEN
                  WRITE(cltxt,*) '   jwrf, jext, jfld: ', jwrf, jext, jfld ,' name: ', TRIM(sndname(jwrf,jext,jfld))
                  CALL wrf_debug(nlevdbg2, cltxt)
               END IF
            END DO
         END DO
      END DO
      CALL wrf_debug(nlevdbg, 'cpl_init: defined variables to be potentially sent' )
      DO jfld = 1, max_cplfld
         DO jext = 1, grid%num_ext_model_couple_dom
            DO jwrf = 1, grid%max_dom
               IF( TRIM(rcvname(jwrf,jext,jfld)) /= 'not defined' ) THEN
                  WRITE(cltxt,*) '   jwrf, jext, jfld: ', jwrf, jext, jfld ,' name: ', TRIM(rcvname(jwrf,jext,jfld))
                  CALL wrf_debug(nlevdbg2, cltxt)
               END IF
            END DO
         END DO
      END DO
      
      ! some checks on grid%cplmask...
      DO jext = 1, grid%num_ext_model_couple_dom

         WRITE(cltxt,*) 'checks on cplmask of external model domain: ', jext               ;   CALL wrf_debug(nlevdbg, cltxt)

         zmin = MINVAL(grid%cplmask(ips:ipe,jext,jps:jpe))
         IF( zmin < 0. ) THEN
            WRITE(cltxt,*) 'min of external model domain cplmask: ',jext,' < 0. : ',zmin   ;   CALL cpl_abort('cpl_defdomain',cltxt)
         END IF
         WRITE(cltxt,*) '   minval(grid%cplmask(ips:ipe,jext,jps:jpe)): ', zmin            ;   CALL wrf_debug(nlevdbg, cltxt)

         zmax = MAXVAL(grid%cplmask(ips:ipe,jext,jps:jpe))
         IF( zmax > 1. ) THEN
            WRITE(cltxt,*) 'max of external model domain cplmask: ',jext,' > 1. : ',zmax   ;   CALL cpl_abort('cpl_defdomain',cltxt)
         END IF
         IF( zmax == 0. ) THEN
            WRITE(cltxt,*) 'max of external model domain cplmask: ',jext,' = 0 '           ;   CALL wrf_message(cltxt)
            WRITE(cltxt,*) '  => no coupling between this external model domain and this WRF patch'   ;   CALL wrf_message(cltxt)
         END IF
         WRITE(cltxt,*) '   maxval(grid%cplmask(ips:ipe,jext,jps:jpe)): ', zmax            ;   CALL wrf_debug(nlevdbg, cltxt)

      END DO
#endif       
      
      IF ( coupler_name == 'oasis' ) CALL cpl_oasis_define( sndname, rcvname, grid )

   END SUBROUTINE cpl_defdomain


   SUBROUTINE cpl_settime( psec )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_settime  ***
      !!
      !! ** Purpose :   update and store the number of second since the beginning of the job.  
      !!--------------------------------------------------------------------
      REAL, INTENT(in) :: psec
      !!--------------------------------------------------------------------

      nsecrun = NINT( psec )
      WRITE(cltxt,*) 'store number of second since the beginning of the job: ', nsecrun   ;   CALL wrf_debug(nlevdbg2, cltxt)

   END SUBROUTINE cpl_settime


   FUNCTION cpl_toreceive( kdomwrf, kdomext, kfldid )
      !!-------------------------------------------------------------------
      !!             ***  FUNCTION cpl_toreceive  ***
      !!
      !! ** Purpose :   send back a logical to tell if a variable must be received or not
      !!--------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER, INTENT(IN) :: kdomext   ! external model domain index
      INTEGER, INTENT(IN) :: kfldid    ! field index
      !
      LOGICAL :: cpl_toreceive
      !!--------------------------------------------------------------------

      IF ( coupler_name == 'oasis' ) cpl_toreceive = cpl_oasis_toreceive( kdomwrf, kdomext, kfldid ) 

   END FUNCTION cpl_toreceive


   FUNCTION cpl_tosend( kdomwrf, kfldid, max_edom )
      !!-------------------------------------------------------------------
      !!             ***  FUNCTION cpl_tosend  ***
      !!
      !! ** Purpose :   send back a logical array to tell if a variable must be
      !!                sent or not to each of the external model domains
      !!--------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER, INTENT(IN) :: kfldid    ! variable index
      INTEGER, INTENT(IN) :: max_edom  ! max number of external model domains
      !
      LOGICAL,DIMENSION(max_edom) :: cpl_tosend
      INTEGER                     :: jext          ! local loop indicees
      !!--------------------------------------------------------------------

      DO jext = 1, max_edom
         IF ( coupler_name == 'oasis' )   cpl_tosend(jext) = cpl_oasis_tosend( kdomwrf, jext, kfldid ) 
      END DO
      
   END FUNCTION cpl_tosend


   FUNCTION cpl_get_fldid( cdsuffix )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_get_fldid  ***
      !!
      !! ** Purpose : send back the field id corresponding to the suffix of a coupling variable name
      !!--------------------------------------------------------------------
      CHARACTER(*), INTENT(IN) :: cdsuffix   ! field name suffix
      !
      INTEGER       :: cpl_get_fldid     ! field index
      INTEGER       :: jfld              ! local loop indicees
      CHARACTER(16) :: clprefix          ! 'WRF_d01_EXT_d01_'
      !!--------------------------------------------------------------------
      cpl_get_fldid = -1   ! default value
         
      clprefix = 'WRF_d01_EXT_d01_'   ! the field id is the same for all WRF domains and external models
      DO jfld = 1, max_cplfld
         IF( clprefix//TRIM(cdsuffix) == sndname(1,1,jfld) )   cpl_get_fldid = jfld
         IF( clprefix//TRIM(cdsuffix) == rcvname(1,1,jfld) )   cpl_get_fldid = jfld
      END DO
          
      IF( cpl_get_fldid == -1 )   CALL cpl_abort( 'cpl_get_fldid', 'variable suffix not found '//TRIM(cdsuffix) )
      WRITE(cltxt,*) 'The id of variable'//TRIM(cdsuffix)//' is: ', cpl_get_fldid   ;   CALL wrf_debug(nlevdbg2, cltxt)

   END FUNCTION cpl_get_fldid

   
   SUBROUTINE cpl_snd( grid )
         !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_snd  ***
      !!
      !! ** Purpose : compute coupling data to be sent and call cpl_sndfield
      !!--------------------------------------------------------------------
      TYPE(domain), INTENT(IN), POINTER :: grid
      !
      INTEGER :: ips,ipe,jps,jpe,kps,kpe ! domain dimension on 1 processor
      INTEGER :: ims,ime,jms,jme,kms,kme ! memory domain dimension on 1 processor 
      INTEGER :: ids,ide,jds,jde,kds,kde ! domain dimension
      INTEGER                           :: ifldid, ifldidu, ifldidv
      INTEGER                           :: i, j
      INTEGER                           :: max_edom
      REAL                              :: zu_uo, zv_vo, zspd, ztmp
      REAL, DIMENSION(:,:), ALLOCATABLE :: cplsnd
      REAL, DIMENSION(:,:), ALLOCATABLE :: taut, taui, tauj
      LOGICAL,DIMENSION(:), ALLOCATABLE :: lltosend
      !!--------------------------------------------------------------------
#if (EM_CORE == 1)
      
      CALL get_ijk_from_grid( grid, ids, ide, jds, jde, kds, kde, &
         &                          ims, ime, jms, jme, kms, kme, &
         &                          ips, ipe, jps, jpe, kps, kpe  )

      max_edom = grid%num_ext_model_couple_dom
      ALLOCATE( cplsnd(ips:ipe,jps:jpe), lltosend(max_edom) )
      
      ! we use ipe and not min(ipe, ide-1) as the variable we are using are coming from grid and are therefore initialized to 0
      !
      !=========================================================================================================
      ! mass fluxes
      !=========================================================================================================
      !
      ! LIQUID_PRECIP : total liquid precipitation (kg/m2/s = mm/s, positive downward)
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'LIQUID_PRECIP' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         cplsnd(ips:ipe,jps:jpe) = ( grid%RAINCV(ips:ipe,jps:jpe) + grid%RAINNCV(ips:ipe,jps:jpe) ) / grid%DT
         CALL cpl_sndfield( grid%id, lltosend, ifldid, cplsnd )
      END IF

      !
      ! SOLID_PRECIP : total solid precipitation (kg/m2/s = mm/s, positive downward)
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'SOLID_PRECIP' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         cplsnd(ips:ipe,jps:jpe) = ( grid%SNOWNCV(ips:ipe,jps:jpe)    &
            &                      + grid%HAILNCV(ips:ipe,jps:jpe) + grid%GRAUPELNCV(ips:ipe,jps:jpe) ) / grid%DT
         CALL cpl_sndfield( grid%id, lltosend, ifldid, cplsnd )
      END IF
      
      ! TOTAL_EVAP : total solid precipitation (kg/m2/s = mm/s, positive upward)
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'TOTAL_EVAP' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         CALL cpl_sndfield( grid%id, lltosend, ifldid, grid%QFX(ips:ipe,jps:jpe) )
      END IF
      
      ! EVAP - PRECIP : total net liquid+solide water fluw (kg/m2/s = mm/s, positive upward)
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'EVAP-PRECIP' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN
         ! first, cumulated liquid+solid precipitation
         cplsnd(ips:ipe,jps:jpe) = grid%RAINCV( ips:ipe,jps:jpe) + grid%RAINNCV(ips:ipe,jps:jpe) &
            &                    + grid%SNOWNCV(ips:ipe,jps:jpe) + grid%HAILNCV(ips:ipe,jps:jpe) + grid%GRAUPELNCV(ips:ipe,jps:jpe)
         ! next, add evaporation
         cplsnd(ips:ipe,jps:jpe) = grid%QFX(ips:ipe,jps:jpe) - cplsnd(ips:ipe,jps:jpe) / grid%DT
         CALL cpl_sndfield( grid%id, lltosend, ifldid, cplsnd )
      END IF
      
      !
      !=========================================================================================================
      ! heat fluxes
      !=========================================================================================================
      !
      ! SURF_NET_SOLAR : total net surface solar heat flux (W.m-2, positive downward)
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'SURF_NET_SOLAR' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         CALL cpl_sndfield( grid%id, lltosend, ifldid, grid%GSW(ips:ipe,jps:jpe) )
      END IF
      
      ! SURF_NET_LONGWAVE : total net surface longwave heat flux (W.m-2, positive downward) : GLW - Thermal_radiation
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'SURF_NET_LONGWAVE' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         cplsnd(ips:ipe,jps:jpe) = grid%GLW(ips:ipe,jps:jpe) - STBOLT * grid%EMISS(ips:ipe,jps:jpe) * grid%SST(ips:ipe,jps:jpe)**4
         CALL cpl_sndfield( grid%id, lltosend, ifldid, cplsnd )
      END IF
      
      ! SURF_LATENT : surface latent heat flux (W.m-2, positive downward) : - LH
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'SURF_LATENT' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         CALL cpl_sndfield( grid%id, lltosend, ifldid, -grid%LH(ips:ipe,jps:jpe) )
      END IF
      
      ! SURF_SENSIBLE : surface sensible heat flux (W.m-2, positive downward) : - HFX
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'SURF_SENSIBLE' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         CALL cpl_sndfield( grid%id, lltosend, ifldid, -grid%HFX(ips:ipe,jps:jpe) )
      END IF
      
      ! SURF_NET_NON-SOLAR : total net surface non-solar heat flux (W.m-2, positive downward)
      !                      GLW - LH - HFX - Thermal_radiation
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'SURF_NET_NON-SOLAR' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         cplsnd(ips:ipe,jps:jpe) = grid%GLW(ips:ipe,jps:jpe) &
            &                      - STBOLT * grid%EMISS(ips:ipe,jps:jpe) * grid%SST(ips:ipe,jps:jpe)**4 &
            &                      - grid%LH(ips:ipe,jps:jpe) - grid%HFX(ips:ipe,jps:jpe) 
         CALL cpl_sndfield( grid%id, lltosend, ifldid, cplsnd )
      END IF
      
      !
      !=========================================================================================================
      ! momentum fluxes
      !=========================================================================================================
      !
      ! TAU[XE] or TAU[YN] or TAUMOD : prepare common fields
      !--------------------------------------------------------------------------------------------------
      ! test if we need to compute the module of the wind speed and stress
      lltosend(:) = cpl_tosend( grid%id, cpl_get_fldid( 'TAUMOD' ), max_edom )   &   ! wind tress module
         &     .OR. cpl_tosend( grid%id, cpl_get_fldid( 'TAUX'   ), max_edom )   &   ! or zonal
         &     .OR. cpl_tosend( grid%id, cpl_get_fldid( 'TAUE'   ), max_edom )       ! or Eastward
      IF ( COUNT(lltosend) > 0 ) THEN 
         ALLOCATE( taut(ips:ipe,jps:jpe), taui(ips:ipe,jps:jpe), tauj(ips:ipe,jps:jpe) )
         DO j = jps, jpe
            DO i = ips, ipe
               zu_uo = grid%u_phy(i,kps,j) - grid%uoce(i,j)
               zv_vo = grid%v_phy(i,kps,j) - grid%voce(i,j)
               zspd = MAX( SQRT( zu_uo**2 + zv_vo**2 ), 1.e-7 )
               taut(i,j) = grid%rho(i,kps,j) * grid%ust(i,j)**2
               ztmp = taut(i,j) / zspd
               taui(i,j) = ztmp * zu_uo
               tauj(i,j) = ztmp * zv_vo
            END DO
         END DO
      END IF
      !
      ! WND_[UE]_01 or WND_[VN]_01 : prepare common fields
      !--------------------------------------------------------------------------------------------------
      lltosend(:) = cpl_tosend( grid%id, cpl_get_fldid( 'WND_U_01' ), max_edom )   &  ! zonal wind speed
         &     .OR. cpl_tosend( grid%id, cpl_get_fldid( 'WND_E_01'   ), max_edom )    ! or eastward
      IF ( COUNT(lltosend) > 0 ) THEN
         DO j = jps, jpe
            DO i = ips, ipe
               zu_uo = grid%u_phy(i,kps,j) - grid%uoce(i,j)
               zv_vo = grid%v_phy(i,kps,j) - grid%voce(i,j)
            END DO
         END DO
      END IF


      ! TAUX / TAUY : Surface zonal / meridional wind stress at atmosphere-ocean interface (N m-2)
      !--------------------------------------------------------------------------------------------------
      ifldidu     = cpl_get_fldid( 'TAUX' )
      ifldidv     = cpl_get_fldid( 'TAUY' )
      lltosend(:) = cpl_tosend( grid%id, ifldidu, max_edom ) .AND. cpl_tosend( grid%id, ifldidv, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         CALL cpl_sndfield( grid%id, lltosend, ifldidu, taui )
         CALL cpl_sndfield( grid%id, lltosend, ifldidv, tauj )
      END IF
      
      ! TAUE / TAUN : Eastward / Northward surface wind stress at atmosphere-ocean interface (N m-2)
      !--------------------------------------------------------------------------------------------------
      ifldidu     = cpl_get_fldid( 'TAUE' )
      ifldidv     = cpl_get_fldid( 'TAUN' )
      lltosend(:) = cpl_tosend( grid%id, ifldidu, max_edom ) .AND. cpl_tosend( grid%id, ifldidv, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         cplsnd(ips:ipe,jps:jpe) = grid%cosa(ips:ipe,jps:jpe) * taui(ips:ipe,jps:jpe)   &
            &                    - grid%sina(ips:ipe,jps:jpe) * tauj(ips:ipe,jps:jpe)
         CALL cpl_sndfield( grid%id, lltosend, ifldidu, cplsnd )
         cplsnd(ips:ipe,jps:jpe) = grid%cosa(ips:ipe,jps:jpe) * tauj(ips:ipe,jps:jpe)   &
            &                    + grid%sina(ips:ipe,jps:jpe) * taui(ips:ipe,jps:jpe)
         CALL cpl_sndfield( grid%id, lltosend, ifldidv, cplsnd )
      END IF
      
      ! TAUMOD : Surface wind stress Module at atmosphere-ocean interface (N m-2)
      !--------------------------------------------------------------------------------------------------
      ifldid      = cpl_get_fldid( 'TAUMOD' )
      lltosend(:) = cpl_tosend( grid%id, ifldid, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN 
         CALL cpl_sndfield( grid%id, lltosend, ifldid, taut )
      END IF

      ! WND_U_01 / WND_V_01 : 1st level zonal / meridional wind speed (m S-1)
      !--------------------------------------------------------------------------------------------------
      ifldidu     = cpl_get_fldid( 'WND_U_01' )
      ifldidv     = cpl_get_fldid( 'WND_V_01' )
      lltosend(:) = cpl_tosend( grid%id, ifldidu, max_edom ) .AND. cpl_tosend( grid%id, ifldidv, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN
         CALL cpl_sndfield( grid%id, lltosend, ifldidu, zu_uo )
         CALL cpl_sndfield( grid%id, lltosend, ifldidv, zv_vo )
      END IF

      ! WND_E_01 / WND_N_01 : Eastward / Northward 1st level wind speed (m s-1)
      !--------------------------------------------------------------------------------------------------
      ifldidu     = cpl_get_fldid( 'WND_E_01' )
      ifldidv     = cpl_get_fldid( 'WND_N_01' )
      lltosend(:) = cpl_tosend( grid%id, ifldidu, max_edom ) .AND. cpl_tosend( grid%id, ifldidv, max_edom )
      IF ( COUNT(lltosend) > 0 ) THEN
         cplsnd(ips:ipe,jps:jpe) = grid%cosa(ips:ipe,jps:jpe) * zu_uo(ips:ipe,jps:jpe)   &
            &                    - grid%sina(ips:ipe,jps:jpe) * zv_vo(ips:ipe,jps:jpe)
         CALL cpl_sndfield( grid%id, lltosend, ifldidu, cplsnd )
         cplsnd(ips:ipe,jps:jpe) = grid%cosa(ips:ipe,jps:jpe) * zv_vo(ips:ipe,jps:jpe)   &
            &                    + grid%sina(ips:ipe,jps:jpe) * zu_uo(ips:ipe,jps:jpe)
         CALL cpl_sndfield( grid%id, lltosend, ifldidv, cplsnd )
      END IF
      
      DEALLOCATE( cplsnd, lltosend )
      IF( ALLOCATED(taut) )   DEALLOCATE(taut, taui,tauj)
#endif       
   END SUBROUTINE cpl_snd


   SUBROUTINE cpl_sndfield( kdomwrf, ldtosend, kfldid, pdata )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_rcv  ***
      !!
      !! ** Purpose :   send coupling data
      !!--------------------------------------------------------------------
      INTEGER,              INTENT(IN) :: kdomwrf   ! wrf domain index
      LOGICAL,DIMENSION(:), INTENT(IN) :: ldtosend
      INTEGER,              INTENT(IN) :: kfldid    ! field index
      REAL, DIMENSION(:,:), INTENT(IN) :: pdata     ! data to be sent
      !
      INTEGER :: jext          ! local loop indicees
      !!--------------------------------------------------------------------

      DO jext = 1, SIZE(ldtosend)
         IF( ldtosend(jext) ) THEN
            IF ( coupler_name == 'oasis' ) CALL cpl_oasis_snd( kdomwrf, jext, kfldid, nsecrun, pdata )
         END IF
      END DO

   END SUBROUTINE cpl_sndfield


   SUBROUTINE cpl_rcv( kdomwrf, cdsuffix,            &
      &                ids, ide, jds, jde, kds, kde, &
      &                ims, ime, jms, jme, kms, kme, &
      &                ips, ipe, jps, jpe, kps, kpe, &
      &                max_edom, pcplmask, pdatacpl, pdataobs )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_rcv  ***
      !!
      !! ** Purpose :   receive coupling data
      !!--------------------------------------------------------------------
      INTEGER,                                                   INTENT(IN   ) :: kdomwrf     ! wrf domain index
      CHARACTER(*),                                              INTENT(IN   ) :: cdsuffix    ! field name suffix
      INTEGER,                                                   INTENT(IN   ) :: ids,ide,jds,jde,kds,kde
      INTEGER,                                                   INTENT(IN   ) :: ims,ime,jms,jme,kms,kme
      INTEGER,                                                   INTENT(IN   ) :: ips,ipe,jps,jpe,kps,kpe
      INTEGER,                                                   INTENT(IN   ) :: max_edom    ! max number of external model domains
      REAL, DIMENSION( ims:ime, 1:max_edom, jms:jme ),           INTENT(IN   ) :: pcplmask    ! coupling mask
      REAL, DIMENSION( ims:ime,             jms:jme ),           INTENT(  OUT) :: pdatacpl    ! coupling data
      REAL, DIMENSION( ims:ime,             jms:jme ), OPTIONAL, INTENT(IN   ) :: pdataobs    ! observed data to be merged
      !
      INTEGER :: jext                                ! external domain index
      INTEGER :: ifldid                              ! field index
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: zdata   ! data received from the coupler
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: obsmsk  ! observation mask
      LOGICAL,DIMENSION(max_edom) :: lltorcv
      !!--------------------------------------------------------------------

      ifldid = cpl_get_fldid( cdsuffix )
      DO jext = 1, max_edom   ! is the variable received from one of the external models?
         lltorcv(jext) = cpl_toreceive( kdomwrf, jext, ifldid )
      END DO

      IF( COUNT(lltorcv) == 0 ) RETURN   ! none of the external models sends this variable
         
      IF( PRESENT(pdataobs) ) THEN
         ! build "observation mask" for this variable
         !    => 1. - SUM( cplmask of involved external model(s) for this specific variable ) 
         obsmsk(ips:ipe,jps:jpe) = 1.0
         DO jext = 1, max_edom
            IF( lltorcv(jext) )   obsmsk(ips:ipe,jps:jpe) = obsmsk(ips:ipe,jps:jpe) - pcplmask(ips:ipe,jext,jps:jpe)
         END DO
         pdatacpl(ips:ipe,jps:jpe) = pdataobs(ips:ipe,jps:jpe) * obsmsk(ips:ipe,jps:jpe)
      ELSE 
         pdatacpl(ips:ipe,jps:jpe) = 0.0   ! default value: consider that pdataobs = 0.
      END IF

      DO jext = 1, max_edom
         IF( lltorcv(jext) ) THEN
            IF( coupler_name == 'oasis' )   CALL cpl_oasis_rcv( kdomwrf, jext, ifldid, nsecrun, zdata )
            pdatacpl(ips:ipe,jps:jpe) = pdatacpl(ips:ipe,jps:jpe) + zdata(ips:ipe,jps:jpe) * pcplmask(ips:ipe,jext,jps:jpe)
         END IF
      END DO

   END SUBROUTINE cpl_rcv
   

   SUBROUTINE cpl_rcv_vector( kdomwrf, cdsuffixu, cdsuffixv, cdsuffixe, cdsuffixn,            &
      &                ids, ide, jds, jde, kds, kde, &
      &                ims, ime, jms, jme, kms, kme, &
      &                ips, ipe, jps, jpe, kps, kpe, &
      &                max_edom, pcplmask, cosa, sina, pdatacplu, pdatacplv, pdataobsu, pdataobsv )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_rcv_vector  ***
      !!
      !! ** Purpose :   receive the 2 components of vector data
      !!                this subroutine is similar to cpl_rcv except that the 2 vector components
      !!                that are received from an external model can be rotated (or not) on the local grid
      !!--------------------------------------------------------------------
      INTEGER,                                                   INTENT(IN   ) :: kdomwrf     ! wrf domain index
      CHARACTER(*),                                              INTENT(IN   ) :: cdsuffixu   ! u-field name suffix
      CHARACTER(*),                                              INTENT(IN   ) :: cdsuffixv   ! v-field name suffix
      CHARACTER(*),                                              INTENT(IN   ) :: cdsuffixe   ! e-field name suffix
      CHARACTER(*),                                              INTENT(IN   ) :: cdsuffixn   ! n-field name suffix
      INTEGER,                                                   INTENT(IN   ) :: ids,ide,jds,jde,kds,kde
      INTEGER,                                                   INTENT(IN   ) :: ims,ime,jms,jme,kms,kme
      INTEGER,                                                   INTENT(IN   ) :: ips,ipe,jps,jpe,kps,kpe
      INTEGER,                                                   INTENT(IN   ) :: max_edom    ! max number of external model domains
      REAL, DIMENSION( ims:ime, 1:max_edom, jms:jme ),           INTENT(IN   ) :: pcplmask    ! coupling mask
      REAL, DIMENSION( ims:ime,             jms:jme ),           INTENT(IN   ) :: cosa        ! Local cosine of map rotation
      REAL, DIMENSION( ims:ime,             jms:jme ),           INTENT(IN   ) :: sina        ! Local   sine of map rotation
      REAL, DIMENSION( ims:ime,             jms:jme ),           INTENT(  OUT) :: pdatacplu   ! coupling u-data on local grid
      REAL, DIMENSION( ims:ime,             jms:jme ),           INTENT(  OUT) :: pdatacplv   ! coupling v-data on local grid
      REAL, DIMENSION( ims:ime,             jms:jme ), OPTIONAL, INTENT(IN   ) :: pdataobsu   ! observed u-data to be merged
      REAL, DIMENSION( ims:ime,             jms:jme ), OPTIONAL, INTENT(IN   ) :: pdataobsv   ! observed v-data to be merged
      !
      INTEGER :: i, j, jext                                   ! loop indexes
      INTEGER :: ifldidu, ifldidv, ifldide, ifldidn           ! field index
      REAL    :: zu, zv                                       ! temporary values
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: zdatau, zdatav   ! data received from the coupler
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: obsmsk           ! observation mask
      LOGICAL,DIMENSION(max_edom) :: lltorcv, lltorcv_uv, lltorcv_en
      !!--------------------------------------------------------------------

      ifldidu = cpl_get_fldid( cdsuffixu )   ;   ifldidv = cpl_get_fldid( cdsuffixv )
      ifldide = cpl_get_fldid( cdsuffixe )   ;   ifldidn = cpl_get_fldid( cdsuffixn )

      DO jext = 1, max_edom   ! is the variable received from one of the external models?
         lltorcv_uv(jext) = cpl_toreceive( kdomwrf, jext, ifldidu ) .AND. cpl_toreceive( kdomwrf, jext, ifldidv )
         lltorcv_en(jext) = cpl_toreceive( kdomwrf, jext, ifldide ) .AND. cpl_toreceive( kdomwrf, jext, ifldidn )
      END DO
      lltorcv = lltorcv_uv .OR. lltorcv_en
      
      IF( COUNT(lltorcv) == 0 ) RETURN   ! none of the external models sends this variable
         
      IF( PRESENT(pdataobsu) .AND. PRESENT(pdataobsv) ) THEN
         ! build "observation mask" for this variable
         !    => 1. - SUM( cplmask of involved EXTERNAL model(s) for this specific variable ) 
         obsmsk(ips:ipe,jps:jpe) = 1.0
         DO jext = 1, max_edom
            IF( lltorcv(jext) )   obsmsk(ips:ipe,jps:jpe) = obsmsk(ips:ipe,jps:jpe) - pcplmask(ips:ipe,jext,jps:jpe)
         END DO
         pdatacplu(ips:ipe,jps:jpe) = pdataobsu(ips:ipe,jps:jpe) * obsmsk(ips:ipe,jps:jpe)
         pdatacplv(ips:ipe,jps:jpe) = pdataobsv(ips:ipe,jps:jpe) * obsmsk(ips:ipe,jps:jpe)
      ELSE 
         pdatacplu(ips:ipe,jps:jpe) = 0.0   ! default value: consider that pdataobsu = 0.
         pdatacplv(ips:ipe,jps:jpe) = 0.0   ! default value: consider that pdataobsv = 0.
      END IF

      DO jext = 1, max_edom
         IF( lltorcv_uv(jext) ) THEN
            IF( coupler_name == 'oasis' ) THEN
               CALL cpl_oasis_rcv( kdomwrf, jext, ifldidu, nsecrun, zdatau )
               CALL cpl_oasis_rcv( kdomwrf, jext, ifldidv, nsecrun, zdatav )
            ENDIF
         ENDIF
         IF( lltorcv_en(jext) ) THEN
            IF( coupler_name == 'oasis' ) THEN
               CALL cpl_oasis_rcv( kdomwrf, jext, ifldide, nsecrun, zdatau )
               CALL cpl_oasis_rcv( kdomwrf, jext, ifldidn, nsecrun, zdatav )
            ENDIF
            DO j = jps, jpe
               DO i = ips, ipe
                  zu = zdatau(i,j)
                  zv = zdatav(i,j)
                  zdatau(i,j) = zu * cosa(i,j) + zv * sina(i,j)
                  zdatav(i,j) = zv * cosa(i,j) - zu * sina(i,j)
               END DO
            END DO
         ENDIF
         IF( lltorcv(jext) ) THEN
            pdatacplu(ips:ipe,jps:jpe) = pdatacplu(ips:ipe,jps:jpe) + zdatau(ips:ipe,jps:jpe) * pcplmask(ips:ipe,jext,jps:jpe)
            pdatacplv(ips:ipe,jps:jpe) = pdatacplv(ips:ipe,jps:jpe) + zdatav(ips:ipe,jps:jpe) * pcplmask(ips:ipe,jext,jps:jpe)
         ENDIF
      END DO

   END SUBROUTINE cpl_rcv_vector


   SUBROUTINE cpl_rcv_sfcdrv( kdomwrf, max_edom, CPLMASK, COSA, SINA,  &
      &                       SST_INPUT, SST, UOCE, VOCE,              &
      &                       CHA_COEF,                                &
      &                       ids, ide, jds, jde, kds, kde,            &
      &                       ims, ime, jms, jme, kms, kme,            &
      &                       ips, ipe, jps, jpe, kps, kpe )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_rcv_sfcdrv  ***
      !!
      !! ** Purpose : reveice fields for the surface driver
      !!
      !-- MAX_EDOM      number of external model domains
      !-- CPLMASK       coupling mask (0 for data read in wrflowinput, 1 data received from the coupler)
      !-- SST           sea-surface temperature (K)
      !-- SST_INPUT     sea-surface temperature read in wrflowinput (K) (= SST if no coupling)
      !-- UOCE          sea surface zonal currents (m s-1)
      !-- VOCE          sea surface meridional currents (m s-1)
      !-- CHA_COEF      charnock coefficient from wave model ()
      !!--------------------------------------------------------------------
      INTEGER,                                         INTENT(IN   ) ::   kdomwrf   ! wrf domain index
      INTEGER,                                         INTENT(IN   ) ::   max_edom    ! max number of external model domains
      INTEGER,                                         INTENT(IN   ) ::   ids,ide,jds,jde,kds,kde
      INTEGER,                                         INTENT(IN   ) ::   ims,ime,jms,jme,kms,kme
      INTEGER,                                         INTENT(IN   ) ::   ips,ipe,jps,jpe,kps,kpe
      REAL, DIMENSION( ims:ime, 1:max_edom, jms:jme ), INTENT(IN   ) ::   CPLMASK
      REAL, DIMENSION( ims:ime,             jms:jme ), INTENT(IN   ) ::   COSA, SINA
      REAL, DIMENSION( ims:ime,             jms:jme ), INTENT(IN   ) ::   SST_INPUT
      REAL, DIMENSION( ims:ime,             jms:jme ), INTENT(  OUT) ::   SST
      REAL, DIMENSION( ims:ime,             jms:jme ), INTENT(  OUT) ::   UOCE
      REAL, DIMENSION( ims:ime,             jms:jme ), INTENT(  OUT) ::   VOCE
      REAL, DIMENSION( ims:ime,             jms:jme ), INTENT(  OUT) ::   CHA_COEF
      !!--------------------------------------------------------------------

      CALL cpl_rcv( kdomwrf, 'SST',                   &
         &              ids, ide, jds, jde, kds, kde, &
         &              ims, ime, jms, jme, kms, kme, &
         &              ips, ipe, jps, jpe, kps, kpe, &
         &              max_edom, CPLMASK, SST, SST_INPUT )
      
      CALL cpl_rcv_vector( kdomwrf, 'UOCE', 'VOCE', 'EOCE', 'NOCE',   &
         &              ids, ide, jds, jde, kds, kde,                 &
         &              ims, ime, jms, jme, kms, kme,                 &
         &              ips, ipe, jps, jpe, kps, kpe,                 &
         &              max_edom, CPLMASK, COSA, SINA, UOCE, VOCE )

      CALL cpl_rcv( kdomwrf, 'CHA_COEF',              &
         &              ids, ide, jds, jde, kds, kde, &
         &              ims, ime, jms, jme, kms, kme, &
         &              ips, ipe, jps, jpe, kps, kpe, &
         &              max_edom, CPLMASK, CHA_COEF )
      
   END SUBROUTINE cpl_rcv_sfcdrv


   SUBROUTINE cpl_store_input( grid, config_flags )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_store_input  ***
      !!
      !! ** Purpose : Store input data that will be merged later with data received from the coupler
      !!--------------------------------------------------------------------
      TYPE(domain)                , INTENT(INOUT) :: grid
      TYPE (grid_config_rec_type) , INTENT(IN   ) :: config_flags
      !
      INTEGER :: ips,ipe,jps,jpe,kps,kpe ! domain dimension on 1 processor
      INTEGER :: ims,ime,jms,jme,kms,kme ! memory domain dimension on 1 processor 
      INTEGER :: ids,ide,jds,jde,kds,kde ! domain dimension
      LOGICAL :: llmust_store
      INTEGER :: jext          ! local loop indicees     
      !!--------------------------------------------------------------------

#if (EM_CORE == 1)
      CALL get_ijk_from_grid( grid, ids, ide, jds, jde, kds, kde, &
         &                          ims, ime, jms, jme, kms, kme, &
         &                          ips, ipe, jps, jpe, kps, kpe  )
      
      ! take care of variables read in AUXINPUT4... 
      ! AUXINPUT4 was just read if:
      ! 1) We asked (legally) for an AUXINPUT4 input AND this is the first time step AFTER an auxinput4_alarm was ringing
      ! OR
      ! 2) This is the first time step
      IF( ( config_flags%auxinput4_interval .NE. 0 .AND. config_flags%io_form_auxinput4 .NE. 0 .AND. grid%just_read_auxinput4 ) &
         .OR. grid%itimestep .EQ. 1 ) THEN
         
         ! if we receive the SST, we need to store the SST that was read in auxinput4 in SST_INPUT
         llmust_store = .FALSE.
         DO jext = 1, grid%num_ext_model_couple_dom
            llmust_store = llmust_store .OR. cpl_toreceive( grid%id, jext, cpl_get_fldid( 'SST' ) )
         END DO
         IF( llmust_store )   grid%sst_input(ips:ipe,jps:jpe) = grid%sst(ips:ipe,jps:jpe)   ! store SST into SST_INPUT 
         
         grid%just_read_auxinput4 = .FALSE.  ! the work as been done and not me done again until we reread data from AUXINPUT4
      
      END IF
#endif      
 
   END SUBROUTINE cpl_store_input


   SUBROUTINE cpl_finalize()
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_finalize  ***
      !!
      !! ** Purpose :   cpl_finalize MPI communications with the coupler
      !!--------------------------------------------------------------------
      IF ( coupler_name == 'oasis' ) CALL cpl_oasis_finalize()

   END SUBROUTINE cpl_finalize


   SUBROUTINE cpl_abort( cdroutine, cdtxt )
      !!-------------------------------------------------------------------
      !!             ***  SUBROUTINE cpl_abort  ***
      !!
      !! ** Purpose :   abort coupling simulation
      !!--------------------------------------------------------------------
      CHARACTER(*), INTENT(IN) :: cdroutine   ! name of the subroutine calling cpl_oasis_abort
      CHARACTER(*), INTENT(IN) :: cdtxt       ! aborting text
      !!--------------------------------------------------------------------

      IF ( coupler_name == 'oasis' ) CALL cpl_oasis_abort( cdroutine, cdtxt )

   END SUBROUTINE cpl_abort


END MODULE module_cpl
