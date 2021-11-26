!WRF:MODEL_LAYER:DYNAMICS
!
MODULE module_damping_em

  USE module_domain             ! frame/module_domain.F
  USE module_wrf_error

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE held_suarez_damp( grid,ru_tend, rv_tend, &
                                rw_tend, t_tend, &
                                rw,t,&
                                ru, rv, p, pb,  &
                               ids,ide, jds,jde, kds,kde, &
                               ims,ime, jms,jme, kms,kme, &
                               its,ite, jts,jte, kts,kte )

    IMPLICIT NONE

   !  Input data.
   TYPE (domain), INTENT(IN) :: grid

    INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
    INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
    INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

    REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),          &
          INTENT(INOUT) ::                         ru_tend, &
                                                   rv_tend, &
                                                   rw_tend, t_tend

    REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),          &
          INTENT(IN) ::                            ru, rv, p, pb, &
                                                           rw,t

    integer :: i,j,k
    real :: time_scale


!     Rayleigh time scale (in seconds)

    call nl_get_asiv_tau_buffer(1,time_scale)

    DO j=jts,min(jde-1,jte)
    DO k=kts,MIN(kte,kde-1)
    DO i=its,ite
         ru_tend(i,k,j)=ru_tend(i,k,j)-grid%buffer_u(i,j)*ru(i,k,j)/time_scale
    END DO
    END DO
    END DO

   DO j = jts,jte
   DO k=kts,MIN(kte,kde-1)
   DO i = its,MIN(ite,ide-1)
          rv_tend(i,k,j)=rv_tend(i,k,j)-grid%buffer_v(i,j)*rv(i,k,j)/time_scale
   END DO
   END DO
   END DO

   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte-1   ! 2,kde   ou kts+1,min(kte,kde-1)  ou max(2,kts),min()
   DO i = its,MIN(ite,ide-1)
          rw_tend(i,k,j) = rw_tend(i,k,j) -grid%buffer_w(i,j)*rw(i,k,j)/time_scale
    END DO
    END DO
    END DO


   DO j = jts,MIN(jte,jde-1)
   DO k = kts,MIN(kte,kde)
   DO i = its,MIN(ite,ide-1)
         t_tend(i,k,j)=t_tend(i,k,j)-grid%buffer_t(i,j)*t(i,k,j)/time_scale
    END DO
    END DO
    END DO


      RETURN
  END SUBROUTINE held_suarez_damp

!------------------------------------------------------------------------------

END MODULE module_damping_em
