subroutine ini_st 
!--------------------------------------------------------------------------------------------------                                              
        USE header
!--------------------------------------------------------------------------------------------------                       
!-------Initial condition with a base stratification with salinity front---------------------------
!--------------------------------------------------------------------------------------------------             
        implicit none 
        integer  i,j,k, wiggles 
        integer, parameter :: npl=832
	REAL(kind=rc_kind) :: tver(npl), sver(npl), dep(npl) 
	REAL(kind=rc_kind) :: bs1(npl), cs1(npl), ds1(npl), bT1(npl), cT1(npl), dT1(npl), &
        & z, seval, g, z1, z2, Tbgrnd, sbgrnd, slfac, slfacnew, beta_s, amplitude, &
        & alpha, thy
                

        open(file="./wig_sf_01/src/z_argo_6901283_new_02.bin", unit=8, access="stream",form="unformatted")
        read(8) dep
        close(8)

        open(file="./wig_sf_01/src/ptemp_argo_6901283_new_02.bin", unit=8, access="stream",form="unformatted")
        read(8) tver
        close(8)

        open(file="./wig_sf_01/src/sal_argo_6901283_new_02.bin", unit=8, access="stream",form="unformatted")
        read(8) sver
        close(8)

      
        mldepth=50d0
        z1=mldepth
        z2=z1+450d0
        !print *, z1, z2		

        wiggles=1
        beta_s=7.4*1d-4
        amplitude=2
        g=9.81d0
        alpha=24d0
        tightness=0.01d0
        slfac=(-1d-3)*alpha/(g*beta_s)
        !print *,slfac
        !yfront=0.5*(yc(NJ)+yc(1))

        !thy=tanh(tightness*(yc(1:NJ)-yfront)*PI)

!       do i=1,npl
!      		print *, dep(i), tver(i), sver(i)
!       end do	
                           
        call spline (npl,dep,sver,bs1,cs1,ds1) 
        call spline (npl,dep,tver,bT1,cT1,dT1) 
                                                                        
        do j=1,NJ 
                do i=1,NI 
                        yfront=0.5d0*(yc(NJ)+yc(1))+amplitude*dsin(2.d0*PI*(xc(i)-xc(1))/(xc(NI)-xc(1))/wiggles)                     
                        thy = tanh(tightness*(yc(j)-yfront)*PI)
                        do k=1,NK
                        z= DL*zc(i,j,k)
                        if (z>-z1) then
                                slfacnew=slfac
                        else if (z.le.-z1 .and. z.ge.-z2) then
                                slfacnew=slfac*(z+z2)**4/(z2-z1)**4
                        else 
                                slfacnew=0
                        end if
                        Tbgrnd = seval(npl,z,dep,tver,bT1,cT1,dT1)
                        sbgrnd = seval(npl,z,dep,sver,bs1,cs1,ds1)
                        T(i,j,k,0)=Tbgrnd
                        s(i,j,k,0)=slfacnew*(1-thy)+sbgrnd      
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
