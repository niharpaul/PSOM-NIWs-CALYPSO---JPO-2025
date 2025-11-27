subroutine ini_st 
!--------------------------------------------------------------------------------------------------                                              
        USE header
!--------------------------------------------------------------------------------------------------                       
!-------Initial condition with a base stratification with density front---------------------------
!--------------------------------------------------------------------------------------------------             
        implicit none 
        integer  i,j,k, wiggles 
        integer, parameter :: npl=803
	REAL(kind=rc_kind) :: sver1(npl), sver2(npl), dep(npl) 
	REAL(kind=rc_kind) :: bs1(npl), cs1(npl), ds1(npl), bs2(npl), cs2(npl), ds2(npl), &
        & z, seval, amplitude, thy, s_south(NK), s_north(NK)
                
        open(file="./cdf_00/src/front_dep.bin", unit=8, access="stream",form="unformatted")
        read(8) dep
        close(8)

        open(file="./cdf_00/src/ss_front_sver.bin", unit=8, access="stream",form="unformatted")
        read(8) sver1
        close(8)
        
        !tver(16:48)=tver(15)

        open(file="./cdf_00/src/ns_front_sver.bin", unit=8,access="stream",form="unformatted")
        read(8) sver2
        close(8)
        
        ! Potential density
        sver1=1000+sver1
        sver2=1000+sver2        
      
        amplitude=1d0
        !print *, z1, z2
	tightness=0.4d0
        !tightness=0.08d0

        wiggles=1
        
        !do i=1,npl
        !print *, dep(i), sver(i)
        !end do	
                           
        call spline (npl,dep,sver1,bs1,cs1,ds1) 
        call spline (npl,dep,sver2,bs2,cs2,ds2)

        do k=1,NK
                z= DL*zc(NI/2,NJ/2,k) 
                s_south(k)= seval(npl,z,dep,sver1,bs1,cs1,ds1)   
                s_north(k)= seval(npl,z,dep,sver2,bs2,cs2,ds2) 
                !write(61,*) k,snorth(k),Tnorth(k)
        end do 
                
        T=0
                                                        
        do j=1,NJ 
                do i=1,NI 
                        yfront=0.5d0*(yc(NJ)+yc(1))+amplitude*dsin(2.d0*PI*(xc(i)-xc(1))/(xc(NI)-xc(1))/wiggles)                     
                        thy = 0.5+ 0.5*tanh(tightness*(yc(j)-yfront)*2.*PI)
                        do k=1,NK
                                s(i,j,k,0)= s_south(k) *(1.d0-thy) + s_north(k)*thy
                        end do 
                end do 
        end do 
!--------------------------------------------------------------------------
!-------Vertical boundary--------------------------------------------------
!--------------------------------------------------------------------------                                                                    
        s(:,:,0,0)=s(:,:,1,0)
        s(:,:,NK+1,0)=s(:,:,NK,0)
        T(:,:,0,0)=T(:,:,1,0)
        T(:,:,NK+1,0)=T(:,:,NK,0)
      
        do k=0,NK+1
!-------------------------------------------------------------------------
!-------NS-boundary-------------------------------------------------------
!-------------------------------------------------------------------------	 
                do i=1,NI 
                        s(i,0,k,0)= s(i,1,k,0) 
                        s(i,NJ+1,k,0)= s(i,NJ,k,0) 
                        T(i,0,k,0)= T(i,1,k,0) 
                        T(i,NJ+1,k,0)= T(i,NJ,k,0) 
                end do 
!-------------------------------------------------------------------------
!-------EW-periodic-------------------------------------------------------
!-------------------------------------------------------------------------                                                        
                do j=0,NJ+1 
                        s(0,j,k,0)= s(NI,j,k,0) 
                        s(NI+1,j,k,0)= s(1,j,k,0) 
                        T(0,j,k,0)= T(NI,j,k,0) 
                        T(NI+1,j,k,0)= T(1,j,k,0) 
                end do 
        end do 
        return 
end subroutine ini_st                                  
