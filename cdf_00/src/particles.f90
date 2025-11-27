#include "cppdefs.h"
#ifdef allow_particle

MODULE particles
  USE header, ONLY : NI,NJ,NK,uf,vf,wf,Jifc,Jjfc,J2d,ux,vy,NPR,wz,PI,dtf,vor,shear,rho,strain,zf,&
                    &s,T,parti_file_num,DL,rc_kind,pcx,pcy,pcz,pcr,dirout,dx,dy,WL,EPS,pfac,&
                    &dztop,D,UL

  ! define the class for particles
  TYPE particle
     REAL(kind=rc_kind) ::  i,j,k,x,y,z,u,v,w,s,t,u0,v0,w0,id,vor,strain,shear,rho,time,divQ,wtotal,wsink,wzf
  END TYPE particle

  TYPE (particle), DIMENSION(:), ALLOCATABLE :: parti
  REAL(kind=rc_kind) ::  dz,swap1,swap2,swap3, mycoef, thedepth, theh, wzf
  INTEGER,ALLOCATABLE :: file_id(:)
  INTEGER :: NPR_eachfile
  CHARACTER(len=3) :: file_id_char 

  PRIVATE   :: NPR_eachfile, file_id_char, dz, swap1, swap2, swap3, file_id, mycoef, thedepth, theh, wzf
  PUBLIC    :: parti

CONTAINS

!==================================================
! OPEN_PARTI_FILES()
!==================================================

! This subroutine divides the total number of particles by the number 
! of particule files specififed in the namelist and creates the necessary files.
  SUBROUTINE open_parti_files()
    IMPLICIT NONE
    INTEGER :: fi

    ALLOCATE (file_id(parti_file_num)) 
    ! Print Error if total number of particles cannot be divided by total number of particle files
    IF (MOD(NPR,parti_file_num) .NE. 0) THEN
      PRINT*, "Error: Please make sure NPR/file_num = integer in mod_particles.f90"
      PRINT*, "Stop model"
      STOP
    ENDIF

    NPR_eachfile = NPR/parti_file_num ! The particle number in each file

    ! Open files
    DO fi = 1, parti_file_num
      WRITE(file_id_char,'(I3.3)') fi
      file_id(fi) = 2000 + fi
      OPEN(file_id(fi), file = TRIM(dirout)//'op.parti-'//file_id_char//'.bin', &
          &form = 'unformatted', access = 'stream', status = 'replace')
    ENDDO

  END SUBROUTINE open_parti_files

!==================================================

!==================================================
! SAVE_PARTI()
!==================================================

! This subroutine writes the listed variables for each of the particle.
  SUBROUTINE save_parti()
    IMPLICIT NONE
    INTEGER :: i_file,ip

    PRINT*,"SAVE PARTICLES"
    SELECT CASE(0)
    CASE (0)
       !save limited variables
       DO i_file = 1, parti_file_num
          DO ip = (i_file - 1) * NPR_eachfile + 1, i_file * NPR_eachfile
             WRITE(file_id(i_file)) parti(ip)%id, &    ! Particle number
                  parti(ip)%time, &          ! time (model timestep)
                  parti(ip)%i, &             ! Position (grid points)
                  parti(ip)%j, &
                  parti(ip)%k, &
                  parti(ip)%x, &             ! Position (distance)
                  parti(ip)%y, &
                  parti(ip)%z, &
                  parti(ip)%u*UL*dx, &             ! Velocity
                  parti(ip)%v*UL*dy, &
                  parti(ip)%w*WL/EPS/parti(ip)%wzf, &             ! Advective velocity [m/s]
                  parti(ip)%wsink*WL/EPS/parti(ip)%wzf, &         ! Sinking velocity [m/s]
                  parti(ip)%wtotal*WL/EPS/parti(ip)%wzf, &        ! Total particle velocity (w + wsink) [m/s]
                  parti(ip)%divQ, &          ! Divergence
                  parti(ip)%rho, &           ! Density
                  parti(ip)%s,   &           ! Salinity
                  parti(ip)%t,   &           ! Temperature
                  parti(ip)%vor, &           ! Vorticity
                  parti(ip)%shear, &         ! Shear
                  parti(ip)%strain           ! Strain
          ENDDO
       ENDDO

    END SELECT

  END SUBROUTINE save_parti
!==================================================


!==================================================
! INI_PARTICLES()
!==================================================

! This subroutine initialize the particle distribution. It first assigns a 0 
! to all variables, then assigns the user-defined distribution to the particles.
  SUBROUTINE ini_particles(time)
    IMPLICIT NONE
    INTEGER :: time,ip,ai

    PRINT*, "initialize files to save particles"
    CALL open_parti_files()

    ! Assign 0-values to all variables
    PRINT*, "# ini particles' velocities",NPR
    DO ip=1, NPR
       parti(ip)%time = DBLE(time)
       parti(ip)%id = DBLE(ip)
       parti(ip)%u0=0d0
       parti(ip)%v0=0d0
       parti(ip)%w0=0d0
       parti(ip)%u=0d0
       parti(ip)%v=0d0
       parti(ip)%w=0d0
       parti(ip)%wsink=0d0
       parti(ip)%wtotal=0d0
       parti(ip)%divQ=0d0
       parti(ip)%s=0d0
       parti(ip)%t=0d0
       parti(ip)%rho=0d0
       parti(ip)%vor=0d0
       parti(ip)%shear=0d0
       parti(ip)%strain=0d0
    ENDDO
    PRINT*, "# finish intial particles' velocities"

    ! User-defined particle positioning.    
    ai = 1.
    DO ip=1, NPR
      parti(ip)%i=REAL(MOD(ip,24)) !REAL(i*REAL(NI)/REAL(NPR))
        IF (MOD(ip,24) .eq. 0) THEN
            parti(ip)%i=24.
            ai = ai+1
        ENDIF
      parti(ip)%j=ai
      parti(ip)%k=20.
    ENDDO
 

  END SUBROUTINE ini_particles
!==================================================


!==================================================
! GET_PARTI_VEL()
!==================================================

! This subroutine interpolates the model fields onto the particles' positions
  SUBROUTINE get_parti_vel(time)
    IMPLICIT NONE
    INTEGER :: i,j,k,ip,ic,jc,kc,ifc,jfc,kfc,time
    REAL(kind=rc_kind) ::  dic,djc,dkc,dif,djf,dkf
    REAL(kind=rc_kind) ::  temp1,temp2,temp3,temp4,vel1,vel2,vel3,vel4,dk1,dk2,dk3,dk4,&
                          &z12,z22,z32,z42,z11,z21,z31,z41,dz1,dz2,dz3,dz4
    REAL(kind=rc_kind), DIMENSION(    0:NI,0:NJ+1        )        :: uxf
    REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ        )        :: vyf
    REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ+1, 0:NK)        :: wzf
    REAL(kind=rc_kind), DIMENSION(      NI,  0:NJ,     NK)        :: vfp
    REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ,   0:NK+1)      :: vf_ex
    REAL(kind=rc_kind), DIMENSION(    0:NI,    NJ,     NK)        :: ufp
    REAL(kind=rc_kind), DIMENSION(    0:NI,  0:NJ+1, 0:NK+1)      :: uf_ex
    REAL(kind=rc_kind), DIMENSION(      NI,    NJ,   0:NK)        :: wfp
    REAL(kind=rc_kind), DIMENSION(    0:NI+1,0:NJ+1, 0:NK)        :: wf_ex
    logical :: exist

    ! COMPUTE WZF
    !== Surface Layer 
    ! Compute wzf in the surface layer by linear interpolation
    DO k = NK-1, NK
        wzf(:,:,k) = 0.5d0*(wz(:,:,k) + wz(:,:,k+1))
    ENDDO
    
    !== Below surface layer
    ! Use the function defined in findz_all
    DO k = 0, NK-2
        DO i = 0, NI+1
            DO j = 0, NJ+1
                mycoef = (D(i,j) +dztop)*exp(pfac)*(exp(-0.5*pfac/(NK-1))-exp(0.5*pfac/(NK-1)))/(exp(pfac)-1)
                wzf(i,j,k) = 1/mycoef*exp(pfac/(NK-1)*parti(ip)%k)
            ENDDO
        ENDDO
    ENDDO

    ! Calculate the face velocities at each k depth-level
    k=0
    wfp(:,:,k) = wf(:,:,k)/J2d(1:NI,1:NJ)*wzf(1:NI,1:NJ,k)
    DO k = 1, NK
       ufp(:,:,k) = uf(:,:,k)/Jifc(:,:,k)
       vfp(:,:,k) = vf(:,:,k)/Jjfc(:,:,k)
       wfp(:,:,k) = wf(:,:,k)/J2d(1:NI,1:NJ)*wzf(1:NI,1:NJ,k)
    ENDDO

    ! Creates extrapolated velocity fields at all three boundaries
    ! Zonal velocity u
    uf_ex=0d0
    uf_ex(:,1:NJ,1:NK) = ufp
    !=== Vertical extrapolation at surface
    uf_ex(:,:,NK+1) = 2*uf_ex(:,:,NK)-uf_ex(:,:,NK-1) ! extrapolation

    ! Meridional velocity v
    vf_ex=0d0
    vf_ex(1:NI,:,1:NK) = vfp
    !=== Zonally periodic
    vf_ex(0,:,:) = vf_ex(NI,:,:)
    vf_ex(NI+1,:,:)=vf_ex(1,:,:)
    !=== Vertical extrapolation at surface
    vf_ex(:,:,NK+1) = 2*vf_ex(:,:,NK)-vf_ex(:,:,NK-1)

    ! Vertical velocity w
    wf_ex=0d0
    wf_ex(1:NI,1:NJ,:) = wfp
    !=== Zonally periodic
    wf_ex(0,:,:) = wf_ex(NI,:,:)
    wf_ex(NI+1,:,:) = wf_ex(1,:,:)
    
    !=== Loops through particles
    DO ip = 1, NPR
       parti(ip)%time=DBLE(time)
       
       IF (parti(ip)%j <= NJ .AND. parti(ip)%j >= 0 .AND. &
            parti(ip)%k <= NK .AND. parti(ip)%k >= 0) THEN
            
          ! =========== WARNING =============
          !ic, jc, kc, is the integer index of the particle relative to 
          !the grids center. Use these values for variables with the ghost points.
          !ifc, jfc, and kfc is the index relative to the coordinates of grid faces. 
          !Use these values for variables on faces.
          
          ic = INT(parti(ip)%i+0.5d0)
          jc = INT(parti(ip)%j+0.5d0)
          kc = INT(parti(ip)%k+0.5d0)

          ifc = INT(parti(ip)%i)
          jfc = INT(parti(ip)%j)
          kfc = INT(parti(ip)%k)

          dif = parti(ip)%i - ifc
          djf = parti(ip)%j - jfc
          dkf = parti(ip)%k - kfc

          dic = parti(ip)%i - ic + 0.5d0
          djc = parti(ip)%j - jc + 0.5d0
          dkc = parti(ip)%k - kc + 0.5d0
          
          !=== Calculate the zonal velocity 
          CALL interp_trilinear(dif,djc,dkc,uf_ex(ifc:ifc+1,jc:jc+1,kc:kc+1),parti(ip)%u)

          !=== Calculate the meridional velocity
          CALL interp_trilinear(dic,djf,dkc,vf_ex(ic:ic+1,jfc:jfc+1,kc:kc+1),parti(ip)%v)

          !=== Calculate the vertical velocity
          CALL interp_trilinear(dic,djc,dkf,wf_ex(ic:ic+1,jc:jc+1,kfc:kfc+1),parti(ip)%w)

          !=== Prescribe a sinking velocity
          ! First, need to determine wzf at the particle's location.
          
          IF (parti(ip)%k<NK-1) THEN
              ! To get wzf, the depth of the water column at the particle's location must be known
              ! The depth of the water column (D) must therefore be interpolated onto the particle's horizontal position.
              ! According to the documentation, D is prescribed at center of i- and j-cells.
              CALL interp_bilinear(dic,djc,D(ic:ic+1,jc:jc+1),thedepth)

              mycoef = (thedepth +dztop)*exp(pfac)*(exp(-0.5*pfac/(NK-1))-exp(0.5*pfac/(NK-1)))/(exp(pfac)-1)
              parti(ip)%wzf = 1/mycoef*exp(pfac/(NK-1)*parti(ip)%k)              
          ELSE
              ! Error if particle released above the NK-1 level because function used 
              ! to determine z-levels is different there and discontinuous with the function used below
              PRINT*, "particles too close to surface layer. Cannot be > NK-1. Code should be fixed."
              stop          
          ENDIF              
          
          ! Then, specify the sinking velocity (in m/s), including the scaling factors
          parti(ip)%wsink= -0d0/86400d0/WL*parti(ip)%wzf*EPS  ! 0 m/day

          ! Finally, compute total vertical velocity
          parti(ip)%wtotal= parti(ip)%w + parti(ip)%wsink

          ! Diagnose other properties
          CALL interp_trilinear(dic,djc,dkc,s(ic:ic+1,jc:jc+1,kc:kc+1,0),parti(ip)%s)
          CALL interp_trilinear(dic,djc,dkc,T(ic:ic+1,jc:jc+1,kc:kc+1,0),parti(ip)%t)
          CALL interp_trilinear(dic,djc,dkc,vor(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%vor)
          CALL interp_trilinear(dic,djc,dkc,rho(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%rho)
          CALL interp_trilinear(dic,djc,dkc,shear(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%shear)
          CALL interp_trilinear(dic,djc,dkc,strain(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%strain)
!          CALL interp_trilinear(dic,djc,dkc,divQ(ic:ic+1,jc:jc+1,kc:kc+1),parti(ip)%divQ)

       ELSE
           !=================      
           ! debug part
           !=================      
           !PRINT*, "particles coordinates are wrong, iPR=",ip,"i,j,k",parti(ip)%i,parti(ip)%j,parti(ip)%k
           !stop
       ENDIF
    ENDDO

  END SUBROUTINE get_parti_vel
!==================================================


!==================================================
! PARTI_FORWARD()
!==================================================

! Based on the particle velocities, the position of the particle is projected in time
! in this subroutine using a 2nd order Adams–Bashforth method
! https://en.wikipedia.org/wiki/Linear_multistep_method

  SUBROUTINE parti_forward()
    IMPLICIT NONE
    INTEGER :: ip

    DO ip = 1, NPR
                
       !=================      
       !=== Assign i-position to particle.
       !=================      
#ifdef periodic_ew
       ! Introduce periodicity in x-direction 
       ! (i.e., particle exiting through east bndry are reseeded in the west, and vice versa)
       
       IF (parti(ip)%i > NI) THEN
          parti(ip)%i = parti(ip)%i - REAL(NI)
       ELSE IF (parti(ip)%i < 0d0) THEN
          parti(ip)%i = parti(ip)%i + REAL(NI)
       ELSE
          parti(ip)%i = parti(ip)%i + 0.5d0 * dtf * (3d0 * parti(ip)%u - parti(ip)%u0)
          
          ! Make sure new position is not outside the domain
          IF (parti(ip)%i > NI) THEN
             parti(ip)%i = parti(ip)%i - REAL(NI)
          ELSE IF (parti(ip)%i < 0d0) THEN
             parti(ip)%i = parti(ip)%i + REAL(NI)
          ENDIF
       ENDIF
       
#else
       ! Damps particle velocity to be 0 at north/south boundaries
       IF (parti(ip)%i > NI-1 .AND. parti(ip)%u > 0) THEN
          parti(ip)%i = parti(ip)%i + parti(ip)%u * dtf / (1d0 + (parti(ip)%u * dtf)/(DBLE(NI) - parti(ip)%i))
       ELSE IF (parti(ip)%i < 1 .AND. parti(ip)%u < 0) THEN
          parti(ip)%i = parti(ip)%i + parti(ip)%u * dtf / (1d0 - dtf/parti(ip)%i)
       ELSE
          parti(ip)%i = parti(ip)%i + 0.5d0 * dtf * (3d0 * parti(ip)%u - parti(ip)%u0)
       ENDIF       
       
#endif
       ! Assign x-position to particle based on model resolution
       parti(ip)%x = parti(ip)%i * dx
       
       !=================      
       !=== Assign j-position to particle.
       !=================      
#ifdef periodic_ns
       ! Introduce periodicity in y-direction 
       ! (i.e., particle exiting through north bndry are reseeded in the south, and vice versa)
       IF (parti(ip)%j > NJ) THEN
          parti(ip)%j = parti(ip)%j - REAL(NJ)
       ELSE IF (parti(ip)%j < 0d0) THEN
          parti(ip)%j = parti(ip)%j + REAL(NJ)
       ELSE
          parti(ip)%j = parti(ip)%j + 0.5d0 * dtf * (3d0 * parti(ip)%v - parti(ip)%v0)
          
          ! Make sure new position is not outside the domain
          IF (parti(ip)%j > NJ) THEN
             parti(ip)%j = parti(ip)%j - REAL(NJ)
          ELSE IF (parti(ip)%j < 0d0) THEN
             parti(ip)%j = parti(ip)%j + REAL(NJ)
          ENDIF
       ENDIF
       
#else
       ! Damps particle velocity to be 0 at north/south boundaries
       IF (parti(ip)%j > NJ-1 .AND. parti(ip)%v > 0) THEN
          parti(ip)%j = parti(ip)%j + parti(ip)%v * dtf / (1d0 + (parti(ip)%v * dtf)/(DBLE(NJ) - parti(ip)%j))
       ELSE IF (parti(ip)%j < 1 .AND. parti(ip)%v < 0) THEN
          parti(ip)%j = parti(ip)%j + parti(ip)%v * dtf / (1d0 - dtf/parti(ip)%j)
       ELSE
          parti(ip)%j = parti(ip)%j + 0.5d0 * dtf * (3d0 * parti(ip)%v - parti(ip)%v0)
       ENDIF
#endif
       ! Assign y-position to particle based on model resolution
       parti(ip)%y = parti(ip)%j * dy
       
       !=================      
       ! Assign k-position to particle
       !=================      
       ! Damps particle velocity to be 0 at surface/bottom boundaries
       IF (parti(ip)%k > NK-1 .AND. parti(ip)%wtotal > 0) THEN
          parti(ip)%k = parti(ip)%k + parti(ip)%wtotal * dtf / (1d0 + (parti(ip)%wtotal * dtf)/(DBLE(NK) - parti(ip)%k))
       !ELSE IF (parti(ip)%k < 1 .AND. parti(ip)%wtotal < 0) THEN
        !  parti(ip)%k = parti(ip)%k + parti(ip)%wtotal * dtf / (1d0 - dtf/parti(ip)%k)
       ELSE
          parti(ip)%k = parti(ip)%k + 0.5d0 * dtf * (3d0 * parti(ip)%wtotal - parti(ip)%w0)
       ENDIF
       
       ! Assign z-position to particle based on sigma level
       ! Calculate the scaled z-depth
       CALL sigma2z(parti(ip)%i,parti(ip)%j,parti(ip)%k,swap1)
       parti(ip)%z = swap1 * DL

      
       !=================      
       !debug part
       !=================      
       IF (parti(ip)%i<0d0 .OR. parti(ip)%i>NI .OR. parti(ip)%j<0d0 &
          .OR. parti(ip)%j>NJ .OR. parti(ip)%k>NK .OR. parti(ip)%k<0d0 ) THEN
          !PRINT*, "particles coordinates are wrong, iPR=",ip,"i,j,k",parti(ip)%i,parti(ip)%j,parti(ip)%k
          !stop
       ENDIF

       !=================      
       ! Replace velocities at previous timestep by new velocities 
       !=================      
       parti(ip)%u0 = parti(ip)%u
       parti(ip)%v0 = parti(ip)%v
       parti(ip)%w0 = parti(ip)%wtotal
    ENDDO

  END SUBROUTINE parti_forward
!==================================================


!==================================================
! INTERP_TRILINEAR()
!==================================================

! This subroutine is used to estimate a variable at a particle's position in 3D
  SUBROUTINE interp_trilinear(di,dj,dk,var,velp)
    !== give 8 corner points of a cube, interpolate point values inside of the cube
    !== di is the distance between the particle to the left face
    !== dj is the distance between the particle to the southern face
    !== dk is the distance between the particle and the bottom face
    IMPLICIT NONE
    REAL(kind=rc_kind), INTENT(in) :: di,dj,dk
    REAL(kind=rc_kind), INTENT(in), DIMENSION(    2,    2  ,   2  )      :: var
    REAL(kind=rc_kind), INTENT(out) :: velp
    REAL(kind=rc_kind) ::  i1,i2,i3,i4,j1,j2

    ! calcuate the Trilinear interpolation
    i1 = (var(2,1,  1)   - var(1,1,  1))*di + var(1,1,  1)
    i2 = (var(2,1,  2) - var(1,1,2))*di + var(1,1,  2)
    i3 = (var(2,2,2) - var(1,2,2))*di +var(1,2,2)
    i4 = (var(2,2,1)   - var(1,2,1))*di + var(1,2,1)

    j1 = (i3 - i2)*dj + i2
    j2 = (i4 - i1)*dj + i1

    velp = (j1 - j2) * dk + j2
  END SUBROUTINE interp_trilinear
  
!==================================================
! INTERP_BILINEAR()
!==================================================

! This subroutine is used to estimate a variable at a particle's position in 2D
  SUBROUTINE interp_bilinear(di,dj,var,velp)
    !== give 4 corner points of a square, interpolate point values inside of the square cell
    !== di is the distance between the particle to the left face
    !== dj is the distance between the particle to the southern face
    IMPLICIT NONE
    REAL(kind=rc_kind), INTENT(in) :: di,dj
    REAL(kind=rc_kind), INTENT(in), DIMENSION(2,2) :: var
    REAL(kind=rc_kind), INTENT(out) :: velp
    REAL(kind=rc_kind) ::  i1,i2

    ! calcuate the Bilinear interpolation
    i1 = (var(2,1)   - var(1,1))*di + var(1,1)
    i2 = (var(2,2)   - var(1,2))*di + var(1,2)

    velp = (i2 - i1)*dj + i1
  END SUBROUTINE interp_bilinear  
    
END MODULE particles

#endif
