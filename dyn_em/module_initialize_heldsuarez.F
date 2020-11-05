!IDEAL:MODEL_LAYER:INITIALIZATION

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!GLmod
!    01/2020 : modif Guillaume Lapeyre pour adapter le vieux code 3.6 en 4.0
!GLend

!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain             ! frame/module_domain.F
   USE module_io_domain          ! share
   USE module_state_description  ! frame
   USE module_model_constants    ! share
   USE module_bc                 ! share
   USE module_timing             ! frame
   USE module_configure          ! frame
   USE module_init_utilities     ! dyn_em
#ifdef DM_PARALLEL
   USE module_dm
#endif

! GLmod ------------------------------------------------------------
  CHARACTER (LEN=256) , PRIVATE :: a_message
! GLend ------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the 
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid 
   !  Local data.
   INTEGER :: idum1, idum2

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include "actual_new_args.inc"
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include "dummy_new_args.inc"
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include "dummy_decl.inc"

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k

  INTEGER :: error
  REAL ::  p_surf, p_level,pd_surf,qvf1, qvf
  REAL :: thtmp, ptmp, temp(3), cof1, cof2,width
  INTEGER :: icm,jcm


! GLmod --------------------------------------------------------------------------

   INTEGER, PARAMETER :: nl_max = 1000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho_local, u, v, qv, pd_in
   INTEGER :: nl_in,im
   INTEGER :: ii, im1, jj, jm1, loop, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, z_at_v, z_at_u
   REAL    :: z_scale, xrad, yrad, zrad, rad, delt
   REAL    :: hm
   REAL    :: pi
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    ::  qvf2_surf,qvf2
   INTEGER :: it,kk

   INTEGER :: stretch_grid
   LOGICAL dry_sounding, debug,with_currents,isbuffer
   INTEGER :: nz_jet, ny_jet, nx_jet
   REAL, DIMENSION(:,:), ALLOCATABLE :: u_jet, rho_jet, th_jet, z_jet
   REAL, DIMENSION(:,:), ALLOCATABLE :: sst_anom,u_o,v_o
   CHARACTER (LEN=256) :: mminlutmp
   REAL, PARAMETER :: htbub=8000., radbub=2000000., radz=8000., tpbub=0.5
   REAL :: piov2, tp,coef_latlon
   INTEGER :: icen, jcen,ig,jg
   REAL :: lbeta,yc,phieq,fcoriolis,beta0
   REAL :: ssteq, sstpole
   REAL :: ysst, lsst
   REAL :: ybeta
   REAL :: arg_phi
   integer :: create_sst
   INTEGER :: is_water
   real :: coef_currents


   real,dimension(:),allocatable :: z_mwr
   integer :: k0
   real :: z_surf,bb_coef,cc_coef,mu_base,z_mean,theta_i,alb_i
   real :: z_top,pression_i
   real :: width_buffer,yy


! GLend   ---------------------------------------------------------------

   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT



! GLmod ------------------------------------------------------------
   piov2 = 2.*atan(1.0)
   icen = ide/4
   jcen = jde/2
   CALL nl_get_asiv_stretch_grid(1,stretch_grid)
   CALL nl_get_asiv_z_scale(1,z_scale)

   delt = 0.

   pi = 2.*asin(1.0)
   nxc = (ide-ids)/4
   nyc = (jde-jds)/2

   nz_jet=kde-kds
   ny_jet=jde-jds
   nx_jet=ide-ids
   ALLOCATE(u_jet(nz_jet,ny_jet))
   ALLOCATE(rho_jet(nz_jet,ny_jet))
   ALLOCATE(th_jet(nz_jet,ny_jet))
   ALLOCATE(z_jet(nz_jet,ny_jet))
   ALLOCATE(sst_anom(nx_jet,ny_jet),u_o(nx_jet,ny_jet),v_o(nx_jet,ny_jet))
! GLend ------------------------------------------------------------


   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

      grid%itimestep=0
      grid%step_number = 0


#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif


! GLmod ------------------------------------------------------------

! est-ce toujours utile de donner des valeurs fausses
!    pour des projections, que nous ne faisons pas ?
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


! Parametrisation de Coriolis
   CALL nl_get_asiv_phieq(1,phieq)
   CALL nl_get_asiv_lbeta(1,lbeta)
   CALL nl_get_asiv_ybeta(1,yc)

   fcoriolis    = 2*eomeg*sin(phieq*pi/180.)
   beta0 = 2*eomeg/6370000*cos(phieq*pi/180.)


! Parametrisation de la SST pour ne pas lire le .bin
   CALL nl_get_asiv_ssteq(1,ssteq)
   CALL nl_get_asiv_sstpole(1,sstpole)
   CALL nl_get_asiv_ysst(1,ysst)
   CALL nl_get_asiv_lsst(1,lsst)
   CALL nl_get_asiv_create_sst(1,create_sst)

   CALL nl_get_asiv_with_currents(1,with_currents)

   CALL nl_get_asiv_isbuffer(1,isbuffer)
   CALL nl_get_asiv_width_buffer(1,width)
   width_buffer=width*1000.
   ! Initialize 2D surface arrays


   coef_latlon=360./2./pi/6370.e3

   DO j = jts, jte
   DO i = its, ite

      ! ig is the I index in the global (domain) span of the array.
      ! jg is the J index in the global (domain) span of the array.
      ig = i - ids + 1  ! ids is not necessarily 1
      jg = j - jds + 1  ! jds is not necessarily 1

! GLmod


!       on centre a phieq=180*y_domain/2/pi/R
      grid%xlat(i,j)  = phieq + &
              coef_latlon*(real(jg)-0.5-real(ide-ids)/2.)*config_flags%dy

! on centre a 180W
      grid%xlong(i,j) = 180. +  &
              coef_latlon*(real(ig)-0.5-real(ide-ids)/2.)*config_flags%dx

      grid%msftx(i,j)    = 1.
      grid%msfux(i,j)    = 1.
      grid%msfvx(i,j)    = 1.
      grid%msfvx_inv(i,j)= 1.

      grid%e(i,j)        = 0.

      if((lbeta>0.).and.(yc>0.).and.(lbeta>0.))then
         grid%f(i,j)        = fcoriolis + beta0*lbeta*1000.*   &
           tanh((config_flags%dy*float(j)-yc*1000)/(lbeta*1000))
      else
         grid%f(i,j)        = fcoriolis
      endif

! GLend --------------------------------------------------


      grid%clat(i,j) = grid%xlat(i,j)
      grid%msfty(i,j) = 1.
      grid%msfuy(i,j) = 1.
      grid%msfvy(i,j) = 1.

      ! The following two are the cosine and sine of the rotation
      ! of projection.  Simple cylindrical is *simple* ... no rotation!
      grid%sina(i,j) = 0.
      grid%cosa(i,j) = 1.

   END DO
   END DO


! GLmod ------------------------------------------------------------



   DO j = jts, jte
   DO i = its, min(ite,ide-1)
      ig = i - ids + 1  
      jg = j - jds + 1  

      grid%xlat_v(i,j)  = phieq + &
              coef_latlon*(real(jg)-1.-real(jde-jds)/2.)*config_flags%dy

      grid%xlong_v(i,j) = 180. + &
              coef_latlon*(real(ig)-0.5-real(ide-ids)/2.)*config_flags%dx

   ENDDO
   ENDDO


   DO j = jts, min(jte,jde-1)
   DO i = its, ite
      ig = i - ids + 1  
      jg = j - jds + 1  

      grid%xlat_u(i,j)  = phieq + &
              coef_latlon*(real(jg)-0.5-real(jde-jds)/2.)*config_flags%dy

      grid%xlong_u(i,j) = 180. +  &
              coef_latlon*(real(ig)-1.-real(ide-ids)/2.)*config_flags%dx

   ENDDO
   ENDDO
! GLend --------------------------------------------------


! ceci va lire la table LANDUSE.TBL
    mminlutmp(1:4)='USGS'
    CALL nl_set_mminlu(1,mminlutmp)
    is_water=16
    grid%iswater=is_water
    CALL nl_set_iswater(1,is_water)


! GLmod --------------------------------------------------

!   DO j=jts,jte
!   DO i=its,ite
!     non staggered
   DO j=jts,min(jte,jde-1)
   DO i=its,min(ite,ide-1)


      ! pour ocean
      grid%xland(i,j)    = 2.  
      grid%landmask(i,j) = 0

      grid%lu_index(i,j) = 16 ! d'apres LANDUSE.TBL
      grid%ivgtyp(i,j)   = 16   ! d'apres VEGPARM.TBL
      grid%isltyp(i,j)   = 14    ! d'apres SOILPARM.TBL

      grid%vegfra(i,j) = 0.
      grid%canwat(i,j) = 0.

      grid%xice(i,j) = 0.
      grid%snow(i,j) = 0.


! ne sert a rien car on lit la table LANDUSE.TBL
!      USGS -> 16eme entree
      grid%albedo(i,j)   = 0.
      grid%emiss(i,j)    = 1.
      grid%thc(i,j)      = 1000. ! thermal inertial

! ne sert a rien non plus
      grid%mavail(i,j)   = 0. 


      grid%znt(i,j)      = 0.01   ! roughness length
                                            ! variera dans le temps...


   END DO
   END DO

! GLend --------------------------------------------------






   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy




!                             profil initial



  write(6,*) ' reading input jet sounding '
  call read_input_jet( u_jet, rho_jet, th_jet, z_jet, nz_jet, ny_jet, &
               grid%ztop, stretch_grid, z_scale )
  write(6,*) ' getting dry sounding for base state '
  write(6,*) ' using middle column in jet sounding, j = ',ny_jet/2
  dry_sounding = .true.

!GLmod
  debug = .false.  !  this will produce print of the sounding
!GLend

  CALL get_sounding( zk, p_in, pd_in, theta, rho_local, u, v, qv, dry_sounding, &
                      nl_max, nl_in, u_jet, rho_jet, th_jet, z_jet,      &
                      nz_jet, ny_jet, ny_jet/2, debug                   )
  write(6,*) ' returned from reading sounding, nl_in is ',nl_in


!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure
!  For the jet, using the middle column for the base state means that
!  we will be extrapolating above the highest height data to the south
!  of the centerline.
 grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )



! GLmod ---------------------------------------------------------------------
    write(6,*) ' --------------------------------------------------------------------- '
    write(6,*) ' pression  au centre du domaine (input_jet) :  ',p_in(1:kte-1)
    write(6,*) ' --------------------------------------------------------------------- '
    write(6,*) ' niveaux z (input_jet) : ',zk(1:kte-1)
    write(6,*) ' --------------------------------------------------------------------- '



!--------------------------------------------------------------------------------
!                 set up the grid

   IF (stretch_grid.eq.0) THEN 

      DO k=1, kde
       grid%znw(k) = 1. - float(k-1)/float(kde-1)
      ENDDO

   ELSEIF((stretch_grid.eq.1).and.(z_scale>0.))THEN

      DO k=1, kde
         grid%znw(k) = - tanh(z_scale*(float(k-1)/float(kde-1)-1.))/tanh(z_scale)
      ENDDO

   ELSEIF((stretch_grid.eq.2).and.(z_scale>0.))THEN

      DO k=1, kde
        grid%znw(k) = (1.-exp(-float(kde-1-(k-1))/float(kde-1)*z_scale))/ &
                    (1.-exp(-z_scale))
      ENDDO


   ELSEIF(stretch_grid.eq.3)THEN


 ! ici on definit les hauteurs du modele
 !     et on s'arrange pour calculer les znw correspondants

      allocate(z_mwr(kde))
      z_mwr=0.

!      premier niveau a 10m, puis espace de 20m <100m
!        avec ~17 niveaux en dessous de 1700m
!         (comme Samelson et al. 2020, MWR)
      k0=4
      do ii=2,k0+4
         z_mwr(ii)=10.+20.*real(ii-2)
      enddo
      z_surf=z_mwr(k0+1)
      

      z_top=grid%ztop
      bb_coef=(20.*real(kde-k0)**2-z_top+z_surf)/real(kde-k0-1)
      cc_coef=z_top-z_surf-bb_coef
      do ii=k0+1,kde-1
         z_mwr(ii+1)=z_surf &
             + bb_coef*real(ii-k0)/real(kde-k0) &
             + cc_coef*real(ii-k0)**2/real(kde-k0)**2
      enddo

      write(6,*) ' '
      write(6,*) ' '
      write(6,*) ' --------------------------------------------------------------------- '
      write(6,*) 'z_mwr : ',z_mwr(1:20)
      write(6,*) 'z_mwr top : ',z_mwr(kde),z_top


! -------------


! idealement on devrait prendre :
    p_surf = interp_0( p_in, zk, 0., nl_in )
    mu_base=p_surf-grid%p_top
!     mais on prend la constante...
!      mu_base=p1000mb-grid%p_top


      grid%znw=1.
      
      do ii=2,kde-1
        
         z_mean=0.5*(z_mwr(ii)+z_mwr(ii-1))

         theta_i=interp_0( theta, zk, z_mean, nl_in )
         pression_i=interp_0( p_in, zk,z_mean, nl_in )
         alb_i=r_d/p1000mb*theta_i*(pression_i/p1000mb)**cvpm

! en fait on inverse la definition des pressions sur les niveaux eta...
        grid%znw(ii)=grid%znw(ii-1) + g*(z_mwr(ii-1)-z_mwr(ii))/mu_base/alb_i

      enddo            

      grid%znw(kde)=0.

      deallocate(z_mwr)


   ELSE
      CALL wrf_message     ( 'ERROR: --- asiv_stretch_grid' )
      CALL wrf_message     ( 'ERROR: --- asiv_stretch_grid=0    ==> Standard linear vertical grid' )
      CALL wrf_message     ( 'ERROR: --- asiv_stretch_grid=1    ==> tanh grid' )
      CALL wrf_message     ( 'ERROR: --- asiv_stretch_grid=2    ==> exp grid' )
      CALL wrf_error_fatal ( 'ERROR: --- Invalid option or z_scale = 0' )
   ENDIF

! GLend ------------------------------------------------------------


   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO

   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)



! --------------------------------------------------
!    ON CALCULE LA GRILLE
! --------------------------------------------------

! Standard WRF Coordinate  : hybrid_opt=0
   DO k=1, kde
      grid%c3f(k) = grid%znw(k)
   ENDDO


   DO k=1, kde
      grid%c4f(k) = ( grid%znw(k) - grid%c3f(k) ) * ( p1000mb - grid%p_top )
   ENDDO


   !  Now on half levels, just add up and divide by 2 (for c3h).  Use (eta-c3)*(p00-pt) for c4 on half levels.
   DO k=1, kde-1
      grid%c3h(k) = ( grid%c3f(k+1) + grid%c3f(k) ) * 0.5
      grid%c4h(k) = ( grid%znu(k) - grid%c3h(k) ) * ( p1000mb - grid%p_top )
   ENDDO



   !  c1 = d(B)/d(eta).  We define c1f as c1 on FULL levels.  For a vertical difference,
   !  we need to use B and eta on half levels.  The k-loop ends up referring to the
   !  full levels, neglecting the top and bottom.

   DO k=kds+1, kde-1
      grid%c1f(k) = ( grid%c3h(k) - grid%c3h(k-1) ) / ( grid%znu(k) - grid%znu(k-1) )
   ENDDO


   !  The boundary conditions to get the coefficients:
   !  1) At k=kts: define d(B)/d(eta) = 1.  This gives us the same value of B and d(B)/d(eta)
   !     when doing the sigma-only B=eta.
   !  2) At k=kte: with the new vertical coordinate, define d(B)/d(eta) = 0.  The curve B SMOOTHLY
   !     goes to zero, and at the very top, B continues to SMOOTHLY go to zero.  Note that for
   !     almost all cases of non B=eta, B is ALREADY=ZERO at the top, so this is a reasonable BC to
   !     assume.
   !  3) At k=kte: when trying to mimic the original vertical coordinate, since B = eta, then
   !     d(B)/d(eta) = 1.

   grid%c1f(kds) = 1.
   grid%c1f(kde) = 1.


   !  c2 = ( 1. - c1(k) ) * (p00 - pt).  There is no vertical differencing, so we can do the
   !  full kds to kde looping.

   DO k=kds, kde
      grid%c2f(k) = ( 1. - grid%c1f(k) ) * ( p1000mb - grid%p_top )
   END DO



   !  Now on half levels for c1 and c2.  The c1h will result from the full level c3 and full
   !  level eta differences.  The c2 value use the half level c1(k).

   DO k=1, kde-1
      grid%c1h(k) = ( grid%c3f(k+1) - grid%c3f(k) ) / ( grid%znw(k+1) - grid%znw(k) )
      grid%c2h(k) = ( 1. - grid%c1h(k) ) * ( p1000mb - grid%p_top )
   END DO



! ---------------------------------------------------------------------
!       on calcule toutes les variables
! ---------------------------------------------------------------------






  DO j=jts,jte
  DO i=its,ite  ! flat surface
    grid%ht(i,j) = 0.
    grid%phb(i,1,j) = 0.  ! g*grid%ht(i,j) en fait
    grid%php(i,1,j) = 0.
    grid%ph0(i,1,j) = 0.  ! grid%phb(i,1,j) en fait
  ENDDO
  ENDDO

! GLend ------------------------------------------------------------

   DO J = jts, jte
   DO I = its, ite

! GLmod ------------------------------------------------------------
      p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
!      p_surf = p0 * EXP(-(g*grid%phb(i,1,j)/g)/(r_d*T0))
! GLend ------------------------------------------------------------

      grid%mub(i,j) = p_surf-grid%p_top

      ! given p (coordinate), calculate theta and compute 1/rho from equation
      ! of state

! GLmod ------------------------------------------------------------

!      DO K = kts, kte-1
    DO K = 1, kte-1

! dans WRF4
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top

!      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top

      grid%pb(i,k,j) = p_level

!         grid%t_init(i,k,j) = T0*(p0/p_level)**rcp
      grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0

      grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
    ENDDO


      ! calculate hydrostatic balance (alternatively we could interpolate
      ! the geopotential from the sounding, but this assures that the base
      ! state is in exact hydrostatic balance with respect to the model eqns.

!      DO k = kts+1, kte
    DO kk  = 2,kte
      k=kk - 1
      grid%phb(i,kk,j) = grid%phb(i,kk-1,j) - grid%dnw(kk-1)*(grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))*grid%alb(i,kk-1,j)
    ENDDO


! subtilite : dans l'ancien code, cela commencait a 2...
!      DO k = kts+1, kte
!         grid%phb(i,k,j) = grid%phb(i,k-1,j) - grid%dnw(k-1)*grid%mub(i,j)*grid%alb(i,k-1,j)
!      ENDDO
! GLend ------------------------------------------------------------

   ENDDO
   ENDDO

   write(6,*) ' --------------------------------------------------------------------- '
   write(6,*) 'phb(1:20)/g : ',grid%phb(1,1:20,1)/g
   write(6,*) 'phb top/g : ',grid%phb(1,kde,1)/g,z_top
   write(6,*) ' --------------------------------------------------------------------- '
   write(6,*) 'znw : ',grid%znw(1:kde-1)
   write(6,*) ' --------------------------------------------------------------------- '
   if(grid%znw(kde-1)<0.)then
     CALL wrf_error_fatal ( 'ERROR: --- some znw < 0.' )
   endif


   ! Now set U & V
   DO J = jts, jte


!   DO I = its, ite
  DO I = its, min(ide-1,ite)

!   DO K = kts, kte-1
  DO K = 1, kte
      grid%u_1(i,k,j) = 0.
      grid%u_2(i,k,j) = 0.
      grid%v_1(i,k,j) = 0.
      grid%v_2(i,k,j) = 0.
   END DO
   END DO
   END DO



   DO j=jts, jte
   DO k=kds, kde
   DO i=its, ite
      grid%ww(i,k,j)  = 0.
      grid%w_1(i,k,j) = 0.
      grid%w_2(i,k,j) = 0.
      grid%h_diabatic(i,k,j) = 0.
   END DO
   END DO
   END DO

! GLmod ----------------------------------------------------------------------
  write(6,*) ' ptop is ',grid%p_top
  write(6,*) ' base state grid%mub(1,1), p_surf is ',grid%mub(1,1),grid%mub(1,1)+grid%p_top
!  calculate full state for each column - this includes moisture.
  write(6,*) ' getting grid%moist sounding for full state '


  dry_sounding = .true.
  IF (config_flags%mp_physics /= 0)  dry_sounding = .false.



   DO im = PARAM_FIRST_SCALAR, num_moist
   DO J = jts, jte
   DO K = kts, kte-1
   DO I = its, ite
      grid%moist(i,k,j,im) = 0.
   END DO
   END DO
   END DO
   END DO


  DO J = jts, min(jde-1,jte)

!  get sounding for this point

  debug = .false.  !  this will turn off print of the sounding
  CALL get_sounding( zk, p_in, pd_in, theta, rho_local, u, v, qv, dry_sounding, &
                      nl_max, nl_in, u_jet, rho_jet, th_jet, z_jet,      &
                      nz_jet, ny_jet, j, debug                          )


 DO I = its, min(ide-1,ite)

!   we could just do the first point in "i" and copy from there, but we'll
!   be lazy and do all the points as if they are all, independent

!   At this point grid%p_top is already set. find the DRY mass in the column 
!   by interpolating the DRY pressure.  

    pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )

!   compute the perturbation mass and the full mass

    grid%mu_1(i,j) = pd_surf-grid%p_top - grid%mub(i,j)
    grid%mu_2(i,j) = grid%mu_1(i,j)
    grid%mu0(i,j) = grid%mu_1(i,j) + grid%mub(i,j)



!   given the dry pressure and coordinate system, interp the potential
!   temperature and qv

! WRF4
!    do k=kds,kde-1
    do k=1,kde-1
!      p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
      p_level = grid%c3h(k)*(pd_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)
    enddo


!   integrate the hydrostatic equation (from the RHS of the bigstep
!   vertical momentum equation) down from the top to get grid%p.
!   first from the top of the model to the top pressure



!WRF4
    kk = kte-1  ! top level
    k=kk+1
    qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2
    grid%p(i,kk,j) = - 0.5*((grid%c1f(k)*grid%Mu_1(i,j))+qvf1*(grid%c1f(k)*grid%Mub(i,j)+grid%c2f(k)))/grid%rdnw(kk)/qvf2
    qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
    grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
    grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)

! seule difference = utilisation des c1f(k)
!    k = kte-1  ! top level
!    qvf1 = 0.5*(grid%moist(i,k,j,P_QV)+grid%moist(i,k,j,P_QV))
!    qvf2 = 1./(1.+qvf1)
!    qvf1 = qvf1*qvf2
!    grid%p(i,k,j) = - 0.5*(grid%mu_1(i,j)+qvf1*grid%mub(i,j))/grid%rdnw(k)/qvf2
!    qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
!    grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
!                (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
!    grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)



!  down the column

!WRF4
    do kk=kte-2,1,-1
      k = kk + 1
      qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%p(i,kk,j) = grid%p(i,kk+1,j) - ((grid%c1f(k)*grid%Mu_1(i,j)) + qvf1*(grid%c1f(k)*grid%Mub(i,j)+grid%c2f(k)))/qvf2/grid%rdn(kk+1)
      qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
      grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                  (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
      grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)
    enddo


! seule difference : utilisation des c1f
!    do k=kte-2,1,-1
!      qvf1 = 0.5*(grid%moist(i,k,j,P_QV)+grid%moist(i,k+1,j,P_QV))
!      qvf2 = 1./(1.+qvf1)
!      qvf1 = qvf1*qvf2
!      grid%p(i,k,j) = grid%p(i,k+1,j) - (grid%mu_1(i,j) + qvf1*grid%mub(i,j))/qvf2/grid%rdn(k+1)
!      qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
!      grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
!                  (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
!      grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
!    enddo

!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.


    grid%ph_1(i,1,j) = 0.

!WRF4
    DO kk  = 2,kte
      k = kk-1
      grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                   ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                    (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )
                                                   
      grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
      grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
    ENDDO


! seule difference, utilisation des c1h
!    dans l'ancien code, commence a 2
!    DO k  = kts+1,kte
!      grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (1./grid%rdnw(k-1))*(       &
!                   (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
!                    grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
!      grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
!      grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
!    ENDDO

! interp u
    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO


  ENDDO
  ENDDO



!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',tpbub

  DO J = jts, min(jde-1,jte)

    yrad = config_flags%dy*float(j-jde/2-1)/radbub

    DO I = its, min(ide-1,ite)

      xrad = float(i-1)/float(ide-ids)

!WRF4
      DO K = 1, kte-1

!     dans l'ancien code, commence a 1
!      DO K = kts, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-htbub)/radz
        RAD=SQRT(yrad*yrad+zrad*zrad)

        IF(RAD <= 1.) THEN

           tp = tpbub*cos(rad*piov2)*cos(rad*piov2)*cos(xrad*2*pi+pi)

           grid%t_1(i,k,j)=grid%t_1(i,k,j)+tp
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically


!WRF4
      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k-1)*grid%mub(i,j)+grid%c2h(k-1))+(grid%c1h(k-1)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k-1)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO


!   write(6,*) ' grid%mu_1 from comp ', grid%mu_1(1,1)
!   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, grid%al, grid%t_1, qv '
!   do k=1,kde-1
!     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1),grid%p(1,k,1), grid%al(1,k,1), &
!                     grid%t_1(1,k,1), grid%moist(1,k,1,P_QV)
!   enddo

!   write(6,*) ' grid%mu_1 from comp ', grid%mu_1(1,1)
!   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
!   write(6,*) ' at j = 1 '
!   do k=1,kde-1
!     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1)+grid%phb(1,k,1), &
!                     grid%p(1,k,1)+grid%pb(1,k,1), grid%al(1,k,1)+grid%alb(1,k,1), &
!                     grid%t_1(1,k,1)+t0, grid%moist(1,k,1,P_QV)
!   enddo


!   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
!   write(6,*) ' at j = jde/2 '
!   do k=1,kde-1
!     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,jde/2)+grid%phb(1,k,jde/2), &
!                     grid%p(1,k,jde/2)+grid%pb(1,k,jde/2), grid%al(1,k,jde/2)+grid%alb(1,k,jde/2), &
!                     grid%t_1(1,k,jde/2)+t0, grid%moist(1,k,jde/2,P_QV)
!   enddo


!   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
!   write(6,*) ' at j = jde-1 '
!   do k=1,kde-1
!     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,jde-1)+grid%phb(1,k,jde-1), &
!                     grid%p(1,k,jde-1)+grid%pb(1,k,jde-1), grid%al(1,k,jde-1)+grid%alb(1,k,jde-1), &
!                     grid%t_1(1,k,jde-1)+t0, grid%moist(1,k,jde-1,P_QV)
!   enddo


!  fill out last i row for u

  DO J = jts, min(jde-1,jte)
  DO I = its, ite

    DO K = 1, kte
      grid%u_1(i,k,j) = grid%u_1(its,k,j)
      grid%u_2(i,k,j) = grid%u_2(its,k,j)
    ENDDO

  ENDDO
  ENDDO


!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO


!    c'est utile un peu partout...
  DO k=1,kte-1
   grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1))/g       

!     ce sont les valeurs pour le Rayleigh damping sur la verticale
    grid%t_base(k) = 0.
    grid%qv_base(k) = 0.
    grid%u_base(k) = 0.
    grid%v_base(k) = 0.
  ENDDO




! Creation du champ de SST
   IF ( create_sst.eq.0 ) THEN
      DO J = jts, min(jde-1,jte)
         DO I = its, min(ide-1,ite)
            grid%tsk(I,J) = (ssteq+sstpole)/2.

           if(with_currents)then
               grid%uoce(I,J)=0.
               grid%voce(I,J)=0.
            endif

         ENDDO
      ENDDO

   ELSEIF(create_sst.eq.1)THEN


       coef_currents=2500*(9.81/1027./1.0e-4)*0.15 &
               *0.5*(ssteq-sstpole)/(lsst*1000.)/2.


      DO J = jts, min(jde-1,jte)
         if(lsst>0.)then
            arg_phi = (config_flags%dy*float(J)/1000.0-ysst)/lsst
         else
            arg_phi=0.
         endif
         if(arg_phi>1.)arg_phi=1.
         if(arg_phi<-1.)arg_phi=-1.
         DO I = its, min(ide-1,ite)

               ! Gradient de SST plus continu entre l equateur et le pole    
               ! telle que dans les modeles de circulation generale
               ! HeldSuarez94 et Frierson05
               grid%tsk(I,J) = (ssteq+sstpole)/2.-(ssteq-sstpole)* &
                    sin(pi*arg_phi/2.)/2.

               if(with_currents)then
                  if((arg_phi>-1.).and.(arg_phi<1.))then
                     grid%uoce(I,J)=coef_currents*cos(pi*arg_phi/2.)
                  else
                      grid%uoce(I,J)=0.
                   endif
                  grid%voce(I,J)=0.
               endif

         ENDDO
      ENDDO

   ELSEIF(create_sst.eq.2)THEN

! on multiplie par un facteur additionnel pour avoir des
!   vorticites de l ordre de 10e-6

      coef_currents=2*5*(2500./pi)*(9.81/1027./1.0e-4)*0.15 &
               *0.5*(ssteq-sstpole)/(lsst*1000.)

      DO J = jts, min(jde-1,jte)
         if(lsst>0.)then 
            arg_phi = (config_flags%dy*float(J)/1000.0-ysst)/lsst
         else
            arg_phi = 0.
         endif
         DO I = its, min(ide-1,ite)

               ! SST cible correspondant a un front plus fort que
               ! B05-B08 pour avoir un front atmosphérique
               grid%tsk(I,J) = ssteq - 0.5*(ssteq - sstpole)*(1+tanh(arg_phi))

               if(with_currents)then
                  grid%uoce(I,J)=coef_currents/cosh(arg_phi)**2
                  grid%voce(I,J)=0.
               endif



         ENDDO
      ENDDO

   ELSEIF (create_sst.eq.3) THEN

      DO J = jts, min(jde-1,jte)
        DO I = its, min(ide-1,ite)
           sst_anom(I,J)=0.
        ENDDO
      ENDDO

      CALL read_sst_field(sst_anom,nx_jet,ny_jet)
   
      DO J = jts, min(jde-1,jte)
       DO I = its, min(ide-1,ite)
        grid%tsk(I,J)=sst_anom(I,J)
       ENDDO
      ENDDO


      if(with_currents)then
         CALL read_velocity_field(u_o,v_o,nx_jet,ny_jet)

         DO J = jts, min(jde-1,jte)
            DO I = its, min(ide-1,ite)
               grid%uoce(I,J)=u_o(I,J)
               grid%voce(I,J)=v_o(I,J)
            ENDDO
         ENDDO
      endif


   ELSE
      CALL wrf_message     ( 'ERROR: --- asiv_create_sst' )
      CALL wrf_message     ( 'ERROR: --- asiv_create_sst=0    ==> constant SST' )
      CALL wrf_message     ( 'ERROR: --- asiv_create_sst=1    ==> Analytic smooth front' )
      CALL wrf_message     ( 'ERROR: --- asiv_create_sst=2    ==> Analytic strong front' )
      CALL wrf_message     ( 'ERROR: --- asiv_create_sst=3    ==> reading from file' )
      CALL wrf_error_fatal ( 'ERROR: --- Invalid option' )
   ENDIF


! tant qu'on y est, on initialise tout a la bonne valeur...
      DO j = jts, MIN(jde-1,jte)
         DO i = its, MIN(ide-1,ite)
            grid%sst(i,j)=grid%tsk(i,j)
            grid%tmn(i,j) = grid%sst(i,j)
         ENDDO
      ENDDO

 DEALLOCATE(u_jet,rho_jet,th_jet,z_jet,sst_anom)


! ------------------------------
!        buffer zone
  if(isbuffer)then

!      width_buffer=jte*config_flags%dy/6.
!      print*,jte*config_flags%dy,width_buffer, (jte*config_flags%dy -width_buffer )

  DO J = jts, min(jde-1,jte)
  DO I = its, ite
      jg = j - jds + 1  ! jds is not necessarily 1

      yy=(real(jg)-0.5)*config_flags%dy

      if(yy<width_buffer)then
         grid%buffer_u(i,J)=1.-yy/width_buffer
      elseif(yy> (jte*config_flags%dy -width_buffer ))then
         grid%buffer_u(i,J)=1. - (jte*config_flags%dy - yy)/width_buffer
      else
         grid%buffer_u(i,J)=0.
      endif

  ENDDO
  ENDDO


  DO J = jts, min(jde-1,jte)
    DO I = its, min(ide-1,ite)
      jg = j - jds + 1  ! jds is not necessarily 1

      yy=(real(jg)-0.5)*config_flags%dy

      if(yy<width_buffer)then
         grid%buffer_t(i,J)=1.-yy/width_buffer
      elseif(yy> (jte*config_flags%dy -width_buffer ))then
         grid%buffer_t(i,J)=1. - (jte*config_flags%dy - yy)/width_buffer
      else
         grid%buffer_t(i,J)=0.
      endif

  ENDDO
  ENDDO



   DO j=jts, jte
   DO i=its, ite

      jg = j - jds + 1          ! jds is not necessarily 1

      yy=(real(jg)-0.5)*config_flags%dy

      if(yy<width_buffer)then
         grid%buffer_w(i,J)=1.-yy/width_buffer
      elseif(yy> (jte*config_flags%dy -width_buffer ))then
         grid%buffer_w(i,J)=1. - (jte*config_flags%dy - yy)/width_buffer
      else
         grid%buffer_w(i,J)=0.
      endif

  ENDDO
  ENDDO

   DO J = jts, jte
   DO I = its, min(ide-1,ite)

      jg = j - jds + 1          ! jds is not necessarily 1

      yy=(real(jg)-1.)*config_flags%dy

      if(yy<width_buffer)then
         grid%buffer_v(i,J)=1.-yy/width_buffer
      elseif(yy> (jte*config_flags%dy -width_buffer ))then
         grid%buffer_v(i,J)=1. - (jte*config_flags%dy - yy)/width_buffer
      else
         grid%buffer_v(i,J)=0.
      endif

  ENDDO
  ENDDO



  endif




      !  Save the dry perturbation potential temperature.
      DO j = jts, min(jde-1,jte)
         DO k = kts, kte
            DO i = its, min(ide-1,ite)
               grid%th_phy_m_t0(i,k,j) = grid%t_2(i,k,j)
            END DO
         END DO
      END DO



!GLend

   RETURN

 END SUBROUTINE init_domain_rk

!---------------------------------------------------------------------

 SUBROUTINE init_module_initialize
 END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

! GL                                --------------------------------------------------
!          ROUTINES SUPPLEMENTAIRES

  subroutine get_sounding( zk, p, p_dry, theta, rho,       &
                             u, v, qv, dry, nl_max, nl_in,  &
                             u_jet, rho_jet, th_jet, z_jet, &
                             nz_jet, ny_jet, j_point, debug )
    implicit none

    integer nl_max, nl_in
    real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
         u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
    logical dry

    integer nz_jet, ny_jet, j_point
    real, dimension(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet

    integer n
    parameter(n=1000)
    logical debug

! input sounding data

    real p_surf, th_surf, qv_surf
    real pi_surf, pi(n)
    real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

    real rho_surf, p_input(n), rho_input(n)
    real pm_input(n)  !  this are for full moist sounding

! local data

    real r
    parameter (r = r_d)
    integer k, it, nl
    real qvf, qvf1, dz

!  first, read the sounding

   call calc_jet_sounding( p_surf, th_surf, qv_surf,                             &
                           h_input, th_input, qv_input, u_input, v_input,        &
                           n, nl, debug, u_jet, rho_jet, th_jet, z_jet, j_point, &
                           nz_jet, ny_jet, dry                                  )

   nl = nz_jet

    if(dry) then
     do k=1,nl
       qv_input(k) = 0.
     enddo
    endif

    if(debug) write(6,*) ' number of input levels = ',nl

      nl_in = nl
      if(nl_in .gt. nl_max ) then
        write(6,*) ' too many levels for input arrays ',nl_in,nl_max
        call wrf_error_fatal ( ' too many levels for input arrays ' )
      end if

    qvf = 1. + rvovrd*qv_input(1) 
    rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
    pi_surf = (p_surf/p1000mb)**(r/cp)

    if(debug) then
      write(6,*) ' surface density is ',rho_surf
      write(6,*) ' surface pi is    ',pi_surf
    end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

        qvf = 1. + rvovrd*qv_input(1) 
        qvf1 = 1. + qv_input(1)
        rho_input(1) = rho_surf
        dz = h_input(1)
        do it=1,10
          pm_input(1) = p_surf &
                  - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
          rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
        enddo

! integrate up the column

        do k=2,nl
          rho_input(k) = rho_input(k-1)
          dz = h_input(k)-h_input(k-1)
          qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
          qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
 
          do it=1,10
            pm_input(k) = pm_input(k-1) &
                    - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
!GLmod
!WRF4
              IF(pm_input(k) .LE. 0. )THEN
                CALL wrf_message("Integrated pressure has gone negative - too cold for chosen height")
                WRITE(a_message,*)'k,pm_input(k),h_input(k),th_input(k) = ',k,pm_input(k),h_input(k),th_input(k)
                CALL wrf_error_fatal ( a_message )
              ENDIF
!GLend

            rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
          enddo
        enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa)   pd(Pa)   theta (K)  den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

     end subroutine get_sounding

!------------------------------------------------------------------

  subroutine calc_jet_sounding( p_surf, th_surf, qv_surf,      &
                                h, th, qv, u, v, n, nl, debug, &
                                u_jet, rho_jet, th_jet, z_jet, &
                                jp, nz_jet, ny_jet, dry       )
  implicit none
  integer :: n, nl, jp, nz_jet, ny_jet

  real, dimension(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet
  real, dimension(n) :: h,th,qv,u,v
  real :: p_surf, th_surf, qv_surf
  logical :: debug, dry

  real, dimension(1:nz_jet) :: rho, rel_hum, p
  integer :: k

!  some local stuff

  real :: tmppi, es, qvs, temperature

!  get sounding from column jp

   do k=1,nz_jet
     h(k)  = z_jet(k,jp)
     th(k) = th_jet(k,jp)
     qv(k) = 0.
     rho(k) = rho_jet(k,jp)
     u(k) = u_jet(k,jp)
     v(k) = 0.
   enddo

   if (.not.dry) then
     DO k=1,nz_jet
       if(h(k) .gt. 8000.) then
         rel_hum(k)=0.1
       else
         rel_hum(k)=(1.-0.90*(h(k)/8000.)**1.25)
       end if
       rel_hum(k) = min(0.7,rel_hum(k))
     ENDDO
   else
     do k=1,nz_jet
       rel_hum(k) = 0.
     enddo
   endif

!  next, compute pressure

   do k=1,nz_jet
     p(k) = p1000mb*(R_d*rho(k)*th(k)/p1000mb)**cpovcv
   enddo

!  here we adjust for fixed moisture profile

     IF (.not.dry)  THEN

!  here we assume the input theta is th_v, so we reset theta accordingly

       DO k=1,nz_jet
         tmppi=(p(k)/p1000mb)**rcp
         temperature = tmppi*th(k)
         if (temperature .gt. svpt0) then
            es  = 1000.*svp1*exp(svp2*(temperature-svpt0)/(temperature-svp3))
            qvs = ep_2*es/(p(k)-es)
         else
            es  = 1000.*svp1*exp( 21.8745584*(temperature-273.16)/(temperature-7.66) )
            qvs = ep_2*es/(p(k)-es)
         endif
         qv(k) = rel_hum(k)*qvs
         th(k) = th(k)/(1.+.61*qv(k))
       ENDDO
 
     ENDIF

!  finally, set the surface data. We'll just do a simple extrapolation

   p_surf = 1.5*p(1) - 0.5*p(2)
   th_surf = 1.5*th(1) - 0.5*th(2)
   qv_surf = 1.5*qv(1) - 0.5*qv(2)

   end subroutine calc_jet_sounding

!---------------------------------------------------------------------


 SUBROUTINE read_input_jet( u, r, t, zk, nz, ny, ztop, stretch_grid, z_scale )

 implicit none

 integer, intent(in) :: nz,ny
 real, dimension(nz,ny), intent(out) :: u,r,t,zk
 real, intent(in) :: ztop,z_scale
 integer, intent(in) :: stretch_grid
 integer :: ny_in, nz_in, j,k
 real, dimension(ny,nz) :: field_in

! this code assumes it is called on processor 0 only

   OPEN(unit=10, file='input_jet', form='unformatted', status='old' )
   REWIND(10) 
   read(10) ny_in,nz_in
   if((ny_in /= ny ) .or. (nz_in /= nz)) then
     write(a_message,*) ' error in input jet dimensions '
     CALL wrf_message (a_message)
     write(a_message,*) ' ny, ny_input, nz, nz_input ', ny, ny_in, nz,nz_in
     CALL wrf_message (a_message)
     write(a_message,*) ' error exit '
     CALL wrf_message (a_message)
     call wrf_error_fatal ( ' error in input jet dimensions ' )
   end if

!GLmod
      write(a_message,*) ' ny, nz ', ny, nz
      CALL wrf_message (a_message)
!GLend

   read(10) field_in
   do j=1,ny
   do k=1,nz
     u(k,j) = field_in(j,k)
   enddo
   enddo

!GLmod
      write(a_message,*) ' max u, min u ', maxval(u),minval(u)
      CALL wrf_message (a_message)
!GLend


   read(10) field_in
   do j=1,ny
   do k=1,nz
     t(k,j) = field_in(j,k)
   enddo
   enddo
!GLmod
      write(a_message,*) ' max t, min t ', maxval(t),minval(t)
      CALL wrf_message (a_message)
!GLend



   read(10) field_in
   do j=1,ny
   do k=1,nz
     r(k,j) = field_in(j,k)
   enddo
   enddo
!GLmod
      write(a_message,*) ' max r, min r ', maxval(r),minval(r)
      CALL wrf_message (a_message)

       write(a_message,*) ' asiv_stretch_grid ', stretch_grid
      CALL wrf_message (a_message)
!GLend



!GLmod
!
!  ATTENTION : il faut que zk soit  en accord avec prepare_init

      do j=1,ny
         DO k=1, nz
            zk(k,j) = ztop*(float(k)-1.)/float(nz-1)
         ENDDO
      enddo



!GLend


 end subroutine read_input_jet


!---------------------------------------------------------------------------

   SUBROUTINE read_sst_field( sst_anom , nx, ny)

   IMPLICIT NONE

   INTEGER, INTENT(IN)                  :: nx,ny
   REAL, DIMENSION(nx,ny), INTENT(OUT)  :: sst_anom
   INTEGER                              :: nx_in, ny_in

   OPEN(UNIT=10, FILE='sst_field.bin', FORM='UNFORMATTED', STATUS='OLD' )
   REWIND(10)
   READ(10) nx_in
   READ(10) ny_in
   IF((nx_in .NE. nx ) .OR. (ny_in .NE. ny)) THEN
     WRITE(0,*) ' error in sst field dimensions '
     WRITE(0,*) ' nx, nx_input, ny, ny_input ', nx, nx_in, ny,ny_in
     WRITE(0,*) ' error exit '
     CALL wrf_error_fatal ( ' error in sst field dimensions ' )
   ENDIF

   READ(10) sst_anom

   END SUBROUTINE read_sst_field


!---------------------------------------------------------------------------

   SUBROUTINE read_velocity_field(u_vel,v_vel , nx, ny)

   IMPLICIT NONE

   INTEGER, INTENT(IN)                  :: nx,ny
   REAL, DIMENSION(nx,ny), INTENT(OUT)  :: u_vel,v_vel
   INTEGER                              :: nx_in, ny_in

   OPEN(UNIT=10, FILE='currents_field.bin', FORM='UNFORMATTED', STATUS='OLD' )
   REWIND(10)
   READ(10) nx_in
   READ(10) ny_in
   IF((nx_in .NE. nx ) .OR. (ny_in .NE. ny)) THEN
     WRITE(0,*) ' error in currents field dimensions '
     WRITE(0,*) ' nx, nx_input, ny, ny_input ', nx, nx_in, ny,ny_in
     WRITE(0,*) ' error exit '
     CALL wrf_error_fatal ( ' error in velocity field dimensions ' )
   ENDIF

   READ(10) u_vel
   READ(10) v_vel

   END SUBROUTINE read_velocity_field





END MODULE module_initialize_ideal
