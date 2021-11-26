MODULE module_xios

   USE module_domain          , ONLY : domain, get_ijk_from_grid
   USE module_configure       , ONLY : grid_config_rec_type
   USE module_driver_constants, ONLY : max_domains
#ifdef key_cpp_xios
   USE xios
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC xios_initmodel, xios_finalizemodel
   PUBLIC xios_initdomain, xios_finalizedomain
   PUBLIC xios_def_field, xios_def_gblatt
   PUBLIC xios_settime, xios_put

#ifdef key_cpp_xios

   LOGICAL, PARAMETER, PUBLIC :: xios_on = .true.

   CHARACTER(256) :: ctxt         ! messages or debug string
   INTEGER, PARAMETER :: nlevdbg0 = 0        ! verbosity level
   INTEGER, PARAMETER :: nlevdbg1 = 10       ! verbosity level
   INTEGER, PARAMETER :: nmaxvar = 1000      ! 
   INTEGER, PARAMETER :: nlenvar = 32      ! 
   INTEGER :: nmpi_comm_wrf       ! MPI communicator between the computing nodes
   INTEGER :: ncxt                ! context id currently used
   INTEGER :: nxiosuniq           ! counter used to define unique ids. 
   LOGICAL, DIMENSION(max_domains) :: lcxtdefopn   ! T if the definition of the context still open
   LOGICAL, PUBLIC   , DIMENSION(max_domains,nmaxvar) :: xios_need
   CHARACTER(nlenvar), DIMENSION(max_domains,nmaxvar) :: cvarnames

   INTERFACE xios_put
      MODULE PROCEDURE xios_p2d, xios_p3d
   END INTERFACE xios_put
   INTERFACE xios_finalizemodel
      MODULE PROCEDURE xios_finalize
   END INTERFACE xios_finalizemodel
   INTERFACE xios_def_gblatt
      MODULE PROCEDURE xios_def_gblatt_char, xios_def_gblatt_int, xios_def_gblatt_real
   END INTERFACE xios_def_gblatt

CONTAINS

   SUBROUTINE xios_initmodel( kl_comm, ldcoupler_on ) 
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE xios_initmodel  ***
      !!
      !! ** Purpose :   initialise WRF-XIOS MPI communications
      !!--------------------------------------------------------------------
      INTEGER, INTENT(INOUT) :: kl_comm       ! local MPI communicator of the model
      LOGICAL, INTENT(IN   ) :: ldcoupler_on  ! .true. if use in coupled mode
      !!--------------------------------------------------------------------
      IF ( ldcoupler_on ) THEN
         CALL xios_initialize( "wrfexe", local_comm = kl_comm )
         WRITE(ctxt,*) 'xios_initmodel: provide mpi comminucator to xios : ', kl_comm
      ELSE
         CALL xios_initialize( "wrfexe", return_comm = kl_comm )
         WRITE(ctxt,*) 'xios_initmodel: get mpi comminucator from xios : ', kl_comm
      ENDIF
      nmpi_comm_wrf = kl_comm
      CALL wrf_debug(nlevdbg1, ctxt)
      ! other default definition
      lcxtdefopn(:) = .TRUE.
      nxiosuniq = 0
      cvarnames(:,:) = "not defined"
      xios_need(:,:) = .FALSE.

   END SUBROUTINE xios_initmodel


   SUBROUTINE xios_initdomain( grid, config_flags )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE xios_initdomain  ***
      !!
      !! ** Purpose :  grid parameters and partitioning 
      !!
      !!----------------------------------------------------------------------
      TYPE(domain), INTENT(IN), POINTER ::   grid
      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
      !
      TYPE(xios_duration) :: dtime    = xios_duration(0, 0, 0, 0, 0, 0)
      TYPE(xios_duration) :: f_op, f_of
      TYPE(xios_date)     :: start_date
      CHARACTER(len=19)   :: cldate
      CHARACTER(len= 7)   :: clwrfdom   ! wrf_d<domain>
      CHARACTER(len=10)   :: tfo
      INTEGER :: ips,ipe,jps,jpe,kps,kpe ! domain dimension on 1 processor
      INTEGER :: ims,ime,jms,jme,kms,kme ! memory domain dimension on 1 processor 
      INTEGER :: ids,ide,jds,jde,kds,kde ! domain dimension
      INTEGER :: iipe, jjpe
      INTEGER :: idn,jdn,ipn,jpn,kpn
      INTEGER :: ji
      INTEGER :: ierr
      !!----------------------------------------------------------------------

      CALL get_ijk_from_grid( grid, ids, ide, jds, jde, kds, kde, &
         &                          ims, ime, jms, jme, kms, kme, &
         &                          ips, ipe, jps, jpe, kps, kpe  )

      WRITE(clwrfdom, fmt="('wrf_d',i2.2)") grid%id
      CALL xios_context_initialize(clwrfdom, nmpi_comm_wrf)
      CALL wrf_debug(nlevdbg1, 'xios_context_initialize '//clwrfdom)
      CALL xios_swap(grid%id)

      ! calendar parameters
      ! +++ get the calendar type ??
      CALL xios_define_calendar( TYPE = "Gregorian", time_origin = xios_date(1900,01,01,00,00,00), &
          &                      start_date = xios_date(config_flags%start_year  , config_flags%start_month ,   &
          &                                             config_flags%start_day   , config_flags%start_hour  ,   &
          &                                             config_flags%start_minute, config_flags%start_second) )
      WRITE(cldate,"(i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)") &
         &   config_flags%start_year, config_flags%start_month , config_flags%start_day   ,    &
         &   config_flags%start_hour, config_flags%start_minute, config_flags%start_second
      CALL wrf_debug(nlevdbg1, 'xios_define_calendar '//cldate//' '//clwrfdom)

      ! horizontal grid definition
      iipe = ipe - COUNT( (/ ipe == ide /) )
      jjpe = jpe - COUNT( (/ jpe == jde /) )
      idn = ide-ids+1   ;   jdn = jde-jds+1

      ipn = iipe-ips+1   ;   jpn = jjpe-jps+1
      CALL set_domain_attr("grid_M", ni_glo = idn-1, nj_glo = jdn-1, ibegin = ips-1, jbegin = jps-1, ni = ipn, nj = jpn)
      CALL set_domain_attr("grid_M", data_dim=2, data_ibegin = 0, data_ni = ipn, data_jbegin = 0, data_nj = jpn)
      CALL set_domain_attr("grid_M", lonvalue = REAL(RESHAPE(grid%XLONG(ips:iipe,jps:jjpe),(/ ipn*jpn /)),kind=8), &
         &                           latvalue = REAL(RESHAPE(grid%XLAT (ips:iipe,jps:jjpe),(/ ipn*jpn /)),kind=8))
      CALL wrf_debug(nlevdbg1, 'grid_M definition ok '//clwrfdom)

      ipn = ipe-ips+1   ;   jpn = jjpe-jps+1
      CALL set_domain_attr("grid_U", ni_glo = idn, nj_glo = jdn-1, ibegin = ips-1, jbegin = jps-1, ni = ipn, nj = jpn)
      CALL set_domain_attr("grid_U", data_dim=2, data_ibegin = 0, data_ni = ipn, data_jbegin = 0, data_nj = jpn)
      CALL set_domain_attr("grid_U", lonvalue = REAL(RESHAPE(grid%XLONG_U(ips:ipe,jps:jjpe),(/ ipn*jpn /)),kind=8),&
         &                           latvalue = REAL(RESHAPE(grid%XLAT_U (ips:ipe,jps:jjpe),(/ ipn*jpn /)),kind=8))
      CALL wrf_debug(nlevdbg1, 'grid_U definition ok '//clwrfdom)

      ipn = iipe-ips+1   ;   jpn = jpe-jps+1
      CALL set_domain_attr("grid_V", ni_glo = idn-1, nj_glo = jdn, ibegin = ips-1, jbegin = jps-1, ni = ipn, nj = jpn)
      CALL set_domain_attr("grid_V", data_dim=2, data_ibegin = 0, data_ni = ipn, data_jbegin = 0, data_nj = jpn)
      CALL set_domain_attr("grid_V", lonvalue = REAL(RESHAPE(grid%XLONG_V(ips:iipe,jps:jpe),(/ ipn*jpn /)),kind=8),&
         &                           latvalue = REAL(RESHAPE(grid%XLAT_V (ips:iipe,jps:jpe),(/ ipn*jpn /)),kind=8))
      CALL wrf_debug(nlevdbg1, 'grid_V definition ok '//clwrfdom)

     ! vertical grid definition
      kpn = kpe-kps+1 - COUNT( (/ kpe == kde /) )
      CALL xios_set_axis_attr("lev_M",n_glo=kpn,VALUE=(/ (REAL(ji,kind=8), ji = 1, kpn) /),   &
         &                            long_name="eta values on half (mass) levels", unit="-", positive="down" )
      CALL wrf_debug(nlevdbg1, 'lev_M definition ok '//clwrfdom)
      kpn = kpe-kps+1
      CALL xios_set_axis_attr("lev_W",n_glo=kpn,VALUE=(/ (REAL(ji,kind=8), ji = 1, kpn) /),   &
         &                            long_name="eta values on full (w) levels", unit="-", positive="down" )
      CALL wrf_debug(nlevdbg1, 'lev_W definition ok '//clwrfdom)

      ! automatic definitions of some of the xml attributs
      ! frequency of the call of xios_put (attribut: freq_op)
      f_op%timestep = 1 ; f_of%timestep = 0
      CALL set_field_attr('field_definition', freq_op=f_op, freq_offset=f_of)
      WRITE(ctxt,*) f_op%timestep, 'ts'
      CALL wrf_debug(nlevdbg1, 'field_definition freq_op '//ctxt)

      ! fields automatic definition: output_wrf will call xios_def_field or each relevant variable
      CALL output_wrf ( 1 , grid , config_flags , 1, ierr )
      ! add attributes to variables used in the outputs files...
      CALL xios_add_field_attr()
      
      ! end file definition
      dtime%second = grid%dt
      CALL xios_set_timestep(dtime)
      CALL wrf_debug(nlevdbg1, 'xios_set_timestep '//clwrfdom)
      CALL xios_close_context_definition()
      lcxtdefopn(grid%id) = .FALSE.
      CALL wrf_debug(nlevdbg1, 'xios_close_context_definition '//clwrfdom)

      CALL xios_update_calendar(0)

      DO ji = 1, nmaxvar
         IF ( TRIM(cvarnames(ncxt,ji)) /= "not defined" )   xios_need(ncxt,ji) = xios_field_is_active( TRIM(cvarnames(ncxt,ji)) )
      END DO
      
   END SUBROUTINE xios_initdomain


   SUBROUTINE xios_finalizedomain( kid )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE xios_finalizedomain  ***
      !!
      !! ** Purpose :  close context associated with this id
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kid       ! domain grid id
      !
      CHARACTER(7) :: clwrfdom   ! wrf_d<domain>
      !!--------------------------------------------------------------------

      WRITE(clwrfdom, fmt="('wrf_d',i2.2)") kid
      IF( xios_is_valid_context(clwrfdom) ) THEN
         CALL xios_swap( kid )
         CALL xios_context_finalize()
         CALL wrf_debug(nlevdbg1, 'xios_finalizedomain '//clwrfdom)
         IF( kid /= 1 ) CALL xios_swap( 1 )   ! return back to wrf_d01 context
      ELSE
         CALL wrf_debug(nlevdbg1, 'xios_finalizedomain not a valid context '//clwrfdom)
      ENDIF

   END SUBROUTINE xios_finalizedomain


   SUBROUTINE xios_swap( kid )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  xios_swap  ***
      !!
      !! ** Purpose :  swap context between different domain grids for xios
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kid       ! domain grid id
      !
      TYPE(xios_context) :: wrf_hdl
      CHARACTER(7)       :: clwrfdom   ! wrf_d<domain>
      !!--------------------------------------------------------------------

      ncxt = kid   ! update current context id
      WRITE(clwrfdom, fmt="('wrf_d',i2.2)") kid
      CALL xios_get_handle(clwrfdom,wrf_hdl)
      CALL xios_set_current_context(wrf_hdl)
      CALL wrf_debug(nlevdbg1, 'xios_swap '//clwrfdom)

   END SUBROUTINE xios_swap


   SUBROUTINE xios_settime( kgrdid, ktime )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE xios_settime   ***
      !!
      !! ** Purpose :  define current timestep calendar
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kgrdid
      INTEGER, INTENT(IN) ::  ktime   
      !!---------------------------------------------------------------------
      IF( kgrdid /= ncxt )   CALL xios_swap(kgrdid)
      CALL xios_update_calendar( ktime )
      WRITE(ctxt,*) 'xios_settime: ', kgrdid, ktime
      CALL wrf_debug(nlevdbg1, ctxt)
      
   END SUBROUTINE xios_settime


   SUBROUTINE xios_def_field( kvarid, cdname, cddesc, cdunit, cdstag, kndim )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE xios_def_field   ***
      !!
      !! ** Purpose :  
      !!---------------------------------------------------------------------
      INTEGER,          INTENT(in) :: kvarid
      CHARACTER(LEN=*), INTENT(in) :: cdname     ! field name
      CHARACTER(LEN=*), INTENT(in) :: cddesc
      CHARACTER(LEN=*), INTENT(in) :: cdunit
      CHARACTER(LEN=*), INTENT(in) :: cdstag
      INTEGER,          INTENT(in) :: kndim
      !
      TYPE(xios_field     ) :: field_hdl
      TYPE(xios_fieldgroup) :: fieldgroup_hdl
      TYPE(xios_variable  ) :: var_hdl
      CHARACTER(128) :: cldesc, clunit
      CHARACTER(  6) :: cldom
      CHARACTER(  5) :: claxe
      CHARACTER( 20) :: cl20
      !!--------------------------------------------------------------------

      IF ( .NOT. xios_is_valid_field(TRIM(cdname)) .AND. lcxtdefopn(ncxt) ) THEN

         ! force default definition
         cl20 = cdname   ! use fixed length to make redeable print
         cldesc = cddesc
         IF( LEN_TRIM(cldesc) == 0 ) cldesc = cdname
         clunit = cdunit
         IF( LEN_TRIM(clunit) == 0 ) clunit = '-'

         cldom = "grid_M"
         IF( SCAN(cdstag,'xX') /= 0 ) cldom = "grid_U"
         IF( SCAN(cdstag,'yY') /= 0 ) cldom = "grid_V"

         CALL xios_get_handle( "field_definition", fieldgroup_hdl)        ! get root group handle
         CALL xios_add_child( fieldgroup_hdl, field_hdl, TRIM(cdname))    ! define new field
         ! define new field attributes
         ! default attributs
         SELECT CASE (kndim)
         CASE(2)
            CALL xios_set_field_attr( cdname, long_name = TRIM(cldesc), unit=TRIM(clunit), domain_ref=cldom )
            CALL wrf_debug(nlevdbg0, 'xios_def_field 2D: '//cl20//' '//TRIM(cldesc))
         CASE(3)
            claxe = "lev_M"
            IF( SCAN(cdstag,'zZ') /=0 ) claxe = "lev_W" 
            CALL xios_set_field_attr( cdname, long_name = TRIM(cldesc), unit=TRIM(clunit), domain_ref=cldom, axis_ref=claxe )
            CALL wrf_debug(nlevdbg0, 'xios_def_field 3D: '//cl20//' '//TRIM(cldesc))
         END SELECT

         IF(           kvarid > nmaxvar )   CALL wrf_error_fatal ( "ERROR: xios_def_field increase nmaxvar" )
         IF( len_TRIM(cdname) > nlenvar )   CALL wrf_error_fatal ( "ERROR: xios_def_field increase nlenvar" )
         cvarnames(ncxt,kvarid) = cdname
         
      ENDIF

   END SUBROUTINE xios_def_field


   SUBROUTINE xios_add_field_attr
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE xios_add_field_attr   ***
      !!
      !! ** Purpose :  
      !!---------------------------------------------------------------------
      TYPE(xios_field)      :: field_hdl
      CHARACTER(len= 14)    :: clfmt      ! write format: ('field',ix.x) 
      CHARACTER(len=  8)    :: clfid      ! field id: fieldx, fieldxx, fieldxxx
      CHARACTER(len=128)    :: clrefid    ! field_ref id (variable name)
      CHARACTER(len=  6)    :: cldom
      CHARACTER(len=  5)    :: claxe
      CHARACTER(len=  3)    :: clstag
      LOGICAL               :: llokref, llokdom, llokaxe
      INTEGER               :: ji, jn                 ! loop counters
      INTEGER               :: jimax                  ! max loop index
      INTEGER               :: imaxdigits = 3         ! 
      !!--------------------------------------------------------------------

      DO jn = 1, imaxdigits
         
         WRITE(clfmt, fmt="('(''field'',i',i1,'.'i1,')')") jn, jn      ! variable atribute format definition            
         jimax = 10**imaxdigits / 10**( imaxdigits - jn ) - 1          ! max loop index
         DO ji = 0, jimax
            
            WRITE(clfid, fmt=clfmt) ji                                 ! field id  
            IF ( xios_is_valid_field(TRIM(clfid)) ) THEN
               CALL xios_is_defined_field_attr( TRIM(clfid), field_ref = llokref)
               IF(llokref) THEN
                 CALL xios_get_field_attr       ( TRIM(  clfid),  field_ref = clrefid )
                 CALL xios_is_defined_field_attr( TRIM(clrefid), domain_ref = llokdom )
                 CALL xios_is_defined_field_attr( TRIM(clrefid),   axis_ref = llokaxe )
                 IF(llokdom) THEN
                    CALL xios_get_field_attr( TRIM(clrefid), domain_ref = cldom )
                    IF(cldom == "grid_M") clstag = ""
                    IF(cldom == "grid_U") clstag = "X"
                    IF(cldom == "grid_V") clstag = "Y"
                    IF(llokaxe) THEN
                       CALL xios_get_field_attr( TRIM(clrefid), axis_ref = claxe )
                       IF(claxe == "lev_W") clstag = TRIM(clstag)//"Z"
                    ENDIF
                    IF( LEN_TRIM(clstag) == 0 ) clstag = '-'
                    CALL set_extra_attr( TRIM(clfid),     "stagger", cvalue = TRIM(clstag) )
                 ENDIF
                 IF(llokaxe) THEN 
                    CALL set_extra_attr( TRIM(clfid), "MemoryOrder", cvalue = "XYZ"        )
                 ELSE
                    CALL set_extra_attr( TRIM(clfid), "MemoryOrder", cvalue = "XY"         )
                 END IF
                    CALL set_extra_attr( TRIM(clfid),   "FieldType", kvalue = 104          )
               ENDIF
            ENDIF
            
         END DO
      END DO
      
   END SUBROUTINE xios_add_field_attr

   !!----------------------------------------------------------------------
   !!                   INTERFACE xios_def_att
   !!----------------------------------------------------------------------

   SUBROUTINE xios_def_gblatt_char( cdname, cvalue )
      CHARACTER(LEN=*), INTENT(in) :: cdname
      CHARACTER(LEN=*), INTENT(in) :: cvalue
      CALL xios_def_gblatt_all( cdname, cvalue = cvalue )
   END SUBROUTINE xios_def_gblatt_char

   SUBROUTINE xios_def_gblatt_int ( cdname, kvalue )
      CHARACTER(LEN=*), INTENT(in) :: cdname
      INTEGER         , INTENT(in) :: kvalue
      CALL xios_def_gblatt_all( cdname, kvalue = kvalue )
   END SUBROUTINE xios_def_gblatt_int

   SUBROUTINE xios_def_gblatt_real( cdname, pvalue )
      CHARACTER(LEN=*), INTENT(in) :: cdname
      REAL            , INTENT(in) :: pvalue
      CALL xios_def_gblatt_all( cdname, pvalue = pvalue )
   END SUBROUTINE xios_def_gblatt_real

   SUBROUTINE xios_def_gblatt_all( cdname, cvalue, kvalue, pvalue )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE xios_def_gblatt_all   ***
      !!
      !! ** Purpose :  
      !!---------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) :: cdname   ! name of this attribute
      CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: cvalue   ! value of this attribute (if character)
      INTEGER         , OPTIONAL, INTENT(in) :: kvalue   ! value of this attribute (if integer)
      REAL            , OPTIONAL, INTENT(in) :: pvalue   ! value of this attribute (if real)
      !
      CHARACTER(len=13) :: clfmt        ! write format: ('file',ix.x)
      CHARACTER(len= 7) :: clfid        ! file id: filex, filexx, filexxx
      INTEGER           :: ji, jn       ! loop counters            
      INTEGER           :: jimax        ! max loop index             
      INTEGER           :: imaxdigits = 3         
      !!--------------------------------------------------------------------

      IF ( lcxtdefopn(ncxt) ) THEN   

         DO jn = 1, imaxdigits
            WRITE(clfmt, fmt="('(''file'',i',i1,'.'i1,')')") jn, jn       ! format definition
            jimax = 10**imaxdigits / 10**( imaxdigits - jn ) - 1          ! max loop index
            DO ji = 0, jimax
               WRITE(clfid, fmt=clfmt) ji                                 ! file id  
               CALL set_extra_attr( TRIM(clfid), cdname, cvalue, kvalue, pvalue, ldglobal = .TRUE.  )
            END DO
         END DO

      ENDIF

   END SUBROUTINE xios_def_gblatt_all


   SUBROUTINE set_extra_attr( cdfid, cdname, cvalue, kvalue, pvalue, ldglobal )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE set_extra_attr   ***
      !!
      !! ** Purpose :  
      !!---------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdfid      ! file/field id
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname     ! name of this attribute  
      CHARACTER(LEN=*), OPTIONAL, INTENT(in) ::   cvalue     ! value of this attribute (if character)
      INTEGER         , OPTIONAL, INTENT(in) ::   kvalue     ! value of this attribute (if integer)  
      REAL            , OPTIONAL, INTENT(in) ::   pvalue     ! value of this attribute (if real)
      LOGICAL         , OPTIONAL, INTENT(in) ::   ldglobal   ! true is global attribute (file attribute instead of field)
      !
      TYPE(xios_file      ) :: file_hdl
      TYPE(xios_field     ) :: field_hdl
      TYPE(xios_variable  ) :: var_hdl
      CHARACTER(LEN=18)     :: clvarid
      LOGICAL               :: llglobal
      LOGICAL               :: llok, llokfile, llokfield
      !!--------------------------------------------------------------------
      llglobal = .FALSE.
      IF(PRESENT(ldglobal))   llglobal = ldglobal

      llokfile  = xios_is_valid_file (TRIM(cdfid)) .AND.       llglobal
      llokfield = xios_is_valid_field(TRIM(cdfid)) .AND. .NOT. llglobal

      IF( llokfile .OR. llokfield ) THEN
         nxiosuniq = nxiosuniq + 1
         WRITE(clvarid, fmt="('unique_xiosid_',i4.4)") nxiosuniq
      ENDIF

      IF(llokfile ) THEN
         CALL xios_get_handle( TRIM(cdfid),  file_hdl          )
         CALL xios_add_child (    file_hdl,   var_hdl, clvarid )
      ENDIF
      IF (llokfield) THEN
         CALL xios_get_handle( TRIM(cdfid), field_hdl          )
         CALL xios_add_child (   field_hdl,   var_hdl, clvarid )
      ENDIF

      IF( llokfile .OR. llokfield ) THEN
         IF( PRESENT(cvalue) ) THEN
            CALL xios_set_attr( var_hdl, TYPE =  "string", name = TRIM(cdname) )
            llok = xios_setVar( clvarid, TRIM(cvalue) )
            WRITE(ctxt,*) cvalue
         ENDIF
         IF( PRESENT(kvalue) ) THEN
            CALL xios_set_attr( var_hdl, TYPE = "int", name = TRIM(cdname) )
            llok = xios_setVar( clvarid,      kvalue  )
            WRITE(ctxt,*) kvalue
         ENDIF
         IF( PRESENT(pvalue) ) THEN
            CALL xios_set_attr( var_hdl, TYPE =   "float", name = TRIM(cdname) )
            llok = xios_setVar( clvarid,      pvalue  )
            WRITE(ctxt,*) pvalue
         ENDIF
         CALL wrf_debug(nlevdbg1, 'set_extra_attr: add attribute '//TRIM(cdname)//' ( => '//TRIM(ctxt)//') to '//TRIM(cdfid))
      ENDIF

   END SUBROUTINE set_extra_attr

   !!----------------------------------------------------------------------
   !!                   INTERFACE xios_put
   !!----------------------------------------------------------------------

   SUBROUTINE xios_p2d( kgrdid, cdname, pfield2d, cdstag   &
      &               , ids,ide,jds,jde,kds,kde    &
      &               , ims,ime,jms,jme,kms,kme    &
      &               , ips,ipe,jps,jpe,kps,kpe    &
      &               , kvarid )
      !!---------------------------------------------------------------------
      INTEGER,                         INTENT(in) :: kgrdid
      CHARACTER(LEN=*),                INTENT(in) :: cdname
      REAL,DIMENSION(ims:ime,jms:jme), INTENT(in) :: pfield2d
      CHARACTER(LEN=*),                INTENT(in) :: cdstag
      INTEGER,                         INTENT(in) :: ids,ide,jds,jde,kds,kde
      INTEGER,                         INTENT(in) :: ims,ime,jms,jme,kms,kme
      INTEGER,                         INTENT(in) :: ips,ipe,jps,jpe,kps,kpe
      INTEGER, OPTIONAL,               INTENT(in) :: kvarid
      !
      INTEGER :: iie,jje
      LOGICAL :: llok
      !!---------------------------------------------------------------------
      IF( kgrdid /= ncxt )   CALL xios_swap(kgrdid)
      IF(PRESENT(kvarid)) THEN
         llok = TRIM(cdname) == TRIM(cvarnames(ncxt,kvarid))
         IF(.NOT.llok) CALL wrf_error_fatal( "xios_p2d: we should not be there! "//   &
            &                                TRIM(cdname)//" /= "//TRIM(cvarnames(ncxt,kvarid)) )
      ELSE
         llok = xios_field_is_active( TRIM(cdname) )
      ENDIF
      
      IF(llok) THEN
         CALL wrf_debug(nlevdbg1, 'xios_p2d: '//TRIM(cdname))
         iie = ipe - COUNT( (/ ipe == ide .AND. SCAN(cdstag,'xX') == 0 /) )
         jje = jpe - COUNT( (/ jpe == jde .AND. SCAN(cdstag,'yY') == 0 /) )
         CALL xios_send_field( cdname, REAL( pfield2d(ips:iie,jps:jje), kind=8 ) )
      ENDIF

   END SUBROUTINE xios_p2d

   SUBROUTINE xios_p3d( kgrdid, cdname, pfield3d, cdstag   &
      &               , ids,ide,jds,jde,kds,kde    &
      &               , ims,ime,jms,jme,kms,kme    &
      &               , ips,ipe,jps,jpe,kps,kpe    &
      &               , kvarid )
      !!---------------------------------------------------------------------
      INTEGER,                                 INTENT(in) :: kgrdid
      CHARACTER(LEN=*),                        INTENT(in) :: cdname
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(in) :: pfield3d
      CHARACTER(LEN=*),                        INTENT(in) :: cdstag
      INTEGER,                                 INTENT(in) :: ids,ide,jds,jde,kds,kde
      INTEGER,                                 INTENT(in) :: ims,ime,jms,jme,kms,kme
      INTEGER,                                 INTENT(in) :: ips,ipe,jps,jpe,kps,kpe
      INTEGER, OPTIONAL,                       INTENT(in) :: kvarid
      !
      INTEGER :: iie,jje,kke
      INTEGER,DIMENSION(3) :: ishape, iorder
      LOGICAL :: llok
      !!---------------------------------------------------------------------
      IF( kgrdid /= ncxt )   CALL xios_swap(kgrdid)
      IF(PRESENT(kvarid)) THEN
         llok = TRIM(cdname) == TRIM(cvarnames(ncxt,kvarid))
         IF(.NOT.llok) CALL wrf_error_fatal( "xios_p3d: we should not be there! "//   &
            &                                TRIM(cdname)//" /= "//TRIM(cvarnames(ncxt,kvarid)) )
      ELSE
         llok = xios_field_is_active( TRIM(cdname) )
      ENDIF

      IF(llok) THEN
         CALL wrf_debug(nlevdbg1, 'xios_p3d: '//TRIM(cdname))
         iie = ipe - COUNT( (/ ipe == ide .AND. SCAN(cdstag,'xX') == 0 /) )
         jje = jpe - COUNT( (/ jpe == jde .AND. SCAN(cdstag,'yY') == 0 /) )
         kke = kpe - COUNT( (/ kpe == kde .AND. SCAN(cdstag,'zZ') == 0 /) )
         ishape = (/ iie - ips + 1, jje - jps + 1, kke - kps + 1 /)
         iorder= (/ 1, 3, 2 /)
         CALL xios_send_field( cdname, REAL( RESHAPE(pfield3d(ips:iie,kps:kke,jps:jje), ishape, order=iorder), kind=8 ) )
      ENDIF

   END SUBROUTINE xios_p3d

   !!----------------------------------------------------------------------
   !!                  set_* routines: interface for groups
   !!----------------------------------------------------------------------

   SUBROUTINE set_domain_attr( cdname, ni_glo, nj_glo, ibegin, jbegin, ni, nj, zoom_ibegin, zoom_jbegin, zoom_ni, zoom_nj,   &
      &                                    data_dim, data_ibegin, data_ni, data_jbegin, data_nj, lonvalue, latvalue )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      INTEGER           , OPTIONAL, INTENT(in) ::   ni_glo, nj_glo, ibegin, jbegin, ni, nj
      INTEGER           , OPTIONAL, INTENT(in) ::   data_dim, data_ibegin, data_ni, data_jbegin, data_nj
      INTEGER           , OPTIONAL, INTENT(in) ::   zoom_ibegin, zoom_jbegin, zoom_ni, zoom_nj
      REAL(kind=8), DIMENSION(:), OPTIONAL, INTENT(in) ::   lonvalue, latvalue

      IF ( xios_is_valid_domain(TRIM(cdname)) ) THEN
         CALL xios_set_domain_attr( cdname, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj  ,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj,   &
            &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, type='curvilinear' )
      ENDIF

      IF ( xios_is_valid_domaingroup(TRIM(cdname)) ) THEN
         CALL xios_set_domaingroup_attr( cdname, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj   ,   &
            &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, type='curvilinear' )
      ENDIF

   END SUBROUTINE set_domain_attr

   SUBROUTINE set_field_attr( cdname, freq_op, freq_offset, long_name, unit )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      TYPE(xios_duration),OPTIONAL , INTENT(in) ::   freq_op
      TYPE(xios_duration),OPTIONAL , INTENT(in) ::   freq_offset
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   long_name
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   unit

      IF ( xios_is_valid_field(TRIM(cdname)) ) THEN
         CALL xios_set_field_attr( cdname, freq_op=freq_op, freq_offset=freq_offset, long_name=long_name, unit=unit )
      ENDIF
      IF ( xios_is_valid_fieldgroup(TRIM(cdname)) ) THEN
         CALL xios_set_fieldgroup_attr( cdname, freq_op=freq_op, freq_offset=freq_offset, long_name=long_name, unit=unit )
      ENDIF

   END SUBROUTINE set_field_attr


   SUBROUTINE set_file_attr( cdname, name, name_suffix )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, name_suffix

      IF ( xios_is_valid_file(TRIM(cdname)) ) THEN
         CALL xios_set_file_attr( cdname, name=name, name_suffix=name_suffix )
      ENDIF
      IF ( xios_is_valid_filegroup(TRIM(cdname)) ) THEN
         CALL xios_set_filegroup_attr( cdname, name=name, name_suffix=name_suffix )
      ENDIF

   END SUBROUTINE set_file_attr


   FUNCTION i2str(int)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: int
      CHARACTER(LEN=255) :: i2str

      WRITE(i2str,*) int
      i2str = ADJUSTL(i2str)

   END FUNCTION i2str

#else

   !!----------------------------------------------------------------------
   !!            DUMMY INTERFACE, not used, needed just to compile...
   !!----------------------------------------------------------------------

   LOGICAL,    PARAMETER, PUBLIC :: xios_on = .FALSE.
   LOGICAL, DIMENSION(1,1), PUBLIC :: xios_need

   INTERFACE xios_def_gblatt
      MODULE PROCEDURE xios_def_gblatt_char, xios_def_gblatt_int, xios_def_gblatt_real
   END INTERFACE xios_def_gblatt
   INTERFACE xios_put
      MODULE PROCEDURE xios_p2d, xios_p3d
   END INTERFACE xios_put

CONTAINS

   SUBROUTINE xios_initmodel( kl_comm, ldcoupler_on ) 
      INTEGER, INTENT(INOUT) :: kl_comm       ! local MPI communicator of the model
      LOGICAL, INTENT(IN   ) :: ldcoupler_on  ! .true. if use in coupled mode
      IF (.FALSE.) WRITE(*,*) kl_comm, ldcoupler_on  ! to avoid compilation warning
   END SUBROUTINE xios_initmodel

   SUBROUTINE xios_initdomain( grid, config_flags )
      TYPE(domain), INTENT(IN), POINTER ::   grid
      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
      IF (.FALSE.) WRITE(*,*) config_flags  ! to avoid compilation warning
   END SUBROUTINE xios_initdomain

   SUBROUTINE xios_finalizedomain( kid )
      INTEGER, INTENT(IN) :: kid       ! domain grid id
      IF (.FALSE.) WRITE(*,*) kid  ! to avoid compilation warning
   END SUBROUTINE xios_finalizedomain

   SUBROUTINE xios_finalizemodel()
      IF (.FALSE.) WRITE(*,*) 'nothing'  ! to avoid compilation warning
   END SUBROUTINE xios_finalizemodel

   SUBROUTINE xios_settime( kgrdid, ktime )
      INTEGER, INTENT(in) ::  kgrdid
      INTEGER, INTENT(in) ::  ktime   
      IF (.FALSE.) WRITE(*,*) kgrdid, ktime  ! to avoid compilation warning      
   END SUBROUTINE xios_settime

   SUBROUTINE xios_def_field( kvarid, cdname, cddesc, cdunit, cdstag, kndim )
      INTEGER,          INTENT(in) :: kvarid
      CHARACTER(LEN=*), INTENT(in) :: cdname
      CHARACTER(LEN=*), INTENT(in) :: cddesc
      CHARACTER(LEN=*), INTENT(in) :: cdunit
      CHARACTER(LEN=*), INTENT(in) :: cdstag
      INTEGER,          INTENT(in) :: kndim
      IF (.FALSE.) WRITE(*,*)  kvarid, cdname, cddesc, cdunit, cdstag, kndim  ! to avoid compilation warning      
   END SUBROUTINE xios_def_field

   SUBROUTINE xios_def_gblatt_char( cdname, cvalue )
      CHARACTER(LEN=*), INTENT(in) :: cdname
      CHARACTER(LEN=*), INTENT(in) :: cvalue
      IF (.FALSE.) WRITE(*,*) cdname, cvalue  ! to avoid compilation warning      
   END SUBROUTINE xios_def_gblatt_char

   SUBROUTINE xios_def_gblatt_int ( cdname, kvalue )
      CHARACTER(LEN=*), INTENT(in) :: cdname
      INTEGER         , INTENT(in) :: kvalue
      IF (.FALSE.) WRITE(*,*) cdname, kvalue  ! to avoid compilation warning      
   END SUBROUTINE xios_def_gblatt_int

   SUBROUTINE xios_def_gblatt_real( cdname, pvalue )
      CHARACTER(LEN=*), INTENT(in) :: cdname
      REAL            , INTENT(in) :: pvalue
      IF (.FALSE.) WRITE(*,*) cdname, pvalue  ! to avoid compilation warning      
   END SUBROUTINE xios_def_gblatt_real

   SUBROUTINE xios_p2d( kgrdid, cdname, pfield2d, cdstag      &
      &               , ids,ide,jds,jde,kds,kde &
      &               , ims,ime,jms,jme,kms,kme &
      &               , ips,ipe,jps,jpe,kps,kpe, kvarid )
      INTEGER,                         INTENT(in) :: kgrdid
      CHARACTER(LEN=*),                INTENT(in) :: cdname
      REAL,DIMENSION(ims:ime,jms:jme), INTENT(in) :: pfield2d
      CHARACTER(LEN=*),                INTENT(in) :: cdstag
      INTEGER,                         INTENT(in) :: ids,ide,jds,jde,kds,kde
      INTEGER,                         INTENT(in) :: ims,ime,jms,jme,kms,kme
      INTEGER,                         INTENT(in) :: ips,ipe,jps,jpe,kps,kpe
      INTEGER, OPTIONAL,               INTENT(in) :: kvarid
      IF (.FALSE.) WRITE(*,*) kgrdid, kvarid, cdname, pfield2d, cdstag  ! to avoid compilation warning      
   END SUBROUTINE xios_p2d

   SUBROUTINE xios_p3d( kgrdid, cdname, pfield3d, cdstag      &
      &               , ids,ide,jds,jde,kds,kde &
      &               , ims,ime,jms,jme,kms,kme &
      &               , ips,ipe,jps,jpe,kps,kpe, kvarid )
      INTEGER,                                 INTENT(in) :: kgrdid
      CHARACTER(LEN=*),                        INTENT(in) :: cdname
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(in) :: pfield3d
      CHARACTER(LEN=*),                        INTENT(in) :: cdstag
      INTEGER,                                 INTENT(in) :: ids,ide,jds,jde,kds,kde
      INTEGER,                                 INTENT(in) :: ims,ime,jms,jme,kms,kme
      INTEGER,                                 INTENT(in) :: ips,ipe,jps,jpe,kps,kpe
      INTEGER, OPTIONAL,                       INTENT(in) :: kvarid
     IF (.FALSE.) WRITE(*,*) kgrdid, kvarid, cdname, pfield3d, cdstag  ! to avoid compilation warning      
   END SUBROUTINE xios_p3d

#endif

END MODULE module_xios
