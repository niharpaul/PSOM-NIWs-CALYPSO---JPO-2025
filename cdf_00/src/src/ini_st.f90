subroutine ini_st 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
!     Across-front use the TANH profile with tightness as a measure of f
!     spread                                                            
!     Larger the factor, tighter the front and larger b_xx.             
!     tightness=4, gives the needed b_xx to make b_xx H/ (2f^2) < 1.    
      implicit none 
      integer  i,j,k,iseed,nplm1,iter,iter2 
      integer, parameter :: npl=832
      REAL(kind=rc_kind) :: rhovert(npl),tvert(npl),svert(npl),dep(npl), &
           depoff(npl) 
      REAL(kind=rc_kind) :: bs1(npl),cs1(npl),ds1(npl),bT1(npl),cT1(npl),dT1(npl), &
      &    z,seval,zmm,sbkgrnd,Tbkgrnd,z1,z2,zprev,alpha,beta1,dtdz,dsdz
      REAL(kind=rc_kind) :: slfac,dscl,rdscl,yarg,ex2y,thy,ran3, &
      &    perturb,slfacnew,dz,bfsqbkgrnd,wiggles,amplitude             
!     tightness = 10 represents a very tight front, =1 loose front      
!     mldepth is the desired ML depth                                   
!     = usual parameter (mldepth= 0.d0, tightness=0.03) ! dy=1km        
!     parameter (mldepth= 100.d0 ) ! tightness below, used in inith_fixd
!     parameter (mldepth= 200.d0 ) ! tightness below, used in inith_fixd
      parameter (alpha=1.6d-4)     ! thermal expansion coeff
      parameter (beta1=7.4d-4)     ! haline coeff	
      	
      open(file="./wig_01/src/z_argo_6901283.bin", unit=8, access="stream",form="unformatted")
      	read(8) dep
      close(8)

      open(file="./wig_01/src/ptemp_argo_6901283.bin", unit=8, access="stream",form="unformatted")
        read(8) tvert
      close(8)

      open(file="./wig_01/src/sal_argo_6901283.bin", unit=8, access="stream",form="unformatted")
        read(8) svert
      close(8)


!      do i=1,832
!       print *, dep(i), tvert(i), svert(i)
!      end do	
                           
!     -----------------                                                 
!     Specify MLD                                                       
      mldepth= 50.d0 
!     mldepth= 200.d0                                                  
                                                                        
!     tightness= 0.3d0    !used in inith_fixdh                          
                          !for larger domain                            
      tightness= 0.03d0 
      bfsqbkgrnd = 1.d-4 
      !dtdz= bfsqbkgrnd/(gpr*10.d0)/alpha
      !dsdz= 0.d0
      dtdz=0
      dsdz=bfsqbkgrnd/(gpr*10.d0)/beta1	

      do k=1,npl 
         depoff(k)= dep(k)- mldepth 
      end do 
                                                                        
      call spline (npl,depoff,svert,bs1,cs1,ds1) 
      call spline (npl,depoff,tvert,bT1,cT1,dT1) 
                                                                        
!     sclwidth is the scale width of the channel (i.e. orig width = 48km
!     yfront is the distance of the front from the coast                
      sclwidth = 48.0 
      yfront = 0.5*(yc(NJ+1) +yc(0)) 
!     -offset yfront = (yc(NJ+1) +yc(0))/3.d0  !used for many MLI runs  
                                                                        
!     ADD a PERTURBATION to the width of the front                      
      iseed= 44294 
      dum = ran3(iseed) 
                                                                        
!     z1 is the depth of the ml (diff rho on both sides above this depth
!     z2 is the vertical extent of the density anamoly (it is gradually 
!        linearly anihillated with depth).                              
!     Orig vals. z1= 50. z2=250.                                        
      z1= mldepth - 50.d0 
!     z2= mldepth + 10.d0 
      z2= mldepth + 50.d0 
!     drho=0.2                      ! for 100-200 m deep ML (using convect)
!     slfac= 0.
      slfac= -(0.3/1025.d0)/beta1 
!     slfac= 0.15d0                 ! for 100 m deep ML (using convect) 
                                                                        
!     0.12 in pot dens, 0.15 in salinity                                
      do j=0,NJ+1 
         do i=0,NI+1 
            do k=0,NK+1 
               z= DL*zc(i,j,k) 
                  if (z.ge.-mldepth) then 
                     Tbkgrnd =                                          &
     &               seval(npl,-1.*mldepth,depoff,tvert,bT1,cT1,dT1)     &
     &                       - (z+mldepth)*dtdz
                     sbkgrnd =                                          &
     &               seval(npl,-1.*mldepth,depoff,svert,bs1,cs1,ds1)     &
     &                       - (z+mldepth)*dsdz
                  else 
                     Tbkgrnd =                                          &
     &                    seval(npl,z,depoff,tvert,bT1,cT1,dT1)        
                     sbkgrnd =                                          &
     &                    seval(npl,z,depoff,svert,bs1,cs1,ds1)        
                  end if 
     
     T(i,j,k,0)=  Tbkgrnd 
     s(i,j,k,0)=  sbkgrnd 

            end do 
         end do 
      end do 
                                                                        
                                                                        
!     WIGGLE in FRONT                                                   
      wiggles=1d0 
      amplitude= 2.0d0 
      do j=0,NJ+1 
                                                                        
         do i=0,NI+1 
            yfront= 0.5d0*(yc(NJ+1) +yc(0)) + amplitude*dsin(2.d0*PI*xc(i)/(0.5*(xc(NI+1)+xc(NI))/wiggles))                     
            thy = tanh(tightness*(yc(j)-yfront)*PI) 
  
            do k=1,NK 
               z= DL*zc(i,j,k) 
                  if (z.ge.-z1) then 
                     slfacnew = slfac 
                  else if (z.ge.-z2) then 
                     slfacnew = slfac*(z+z2)/(z2-z1) 
                  else 
                     slfacnew = 0.d0 
                  end if 
!                  if ((i.eq.1).and.(j.eq.1)) write(6,*) 'slfc',k,slfacnew
                                                                        
                  !T(i,j,k,0)= slfacnew*(thy-1.d0) + T(i,NJ,k,0) 
                   s(i,j,k,0)= slfacnew*(-thy+1.d0) + s(i,NJ,k,0)                                                     
            end do 
         end do 
      end do 
!      write(6,*) 'rho i=24,k=24', (s(24,j,24,0),j=1,NJ)                
!      stop                                                             
                                                                        
!         do iter2=1,10                                                 
!            call conadjust(0)                                          
!         end do                                                        
      
      do k=0,NK+1 
         do i=1,NI 
            s(i,0,k,0)= s(i,1,k,0) 
            s(i,NJ+1,k,0)= s(i,NJ,k,0) 
            T(i,0,k,0)= T(i,1,k,0) 
            T(i,NJ+1,k,0)= T(i,NJ,k,0) 
      end do 
!     periodicew                                                        
         do j=0,NJ+1 
            s(0,j,k,0)= s(NI,j,k,0) 
            s(NI+1,j,k,0)= s(1,j,k,0) 
            T(0,j,k,0)= T(NI,j,k,0) 
            T(NI+1,j,k,0)= T(1,j,k,0) 
         end do 
      end do 
                                                                        
                                                                        
      return 
end subroutine ini_st                                        
