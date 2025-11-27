subroutine wind_stress(udif,vdif,step) 
  !     ---------------------------------------------                     
  USE header

  implicit none 
  integer i,j,k,m,step
                                                                        
  REAL(kind=rc_kind) :: udif(NI,NJ,NK), vdif(NI,NJ,NK), stressprofile(NJ) 

  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv, amplitude

  REAL(kind=rc_kind)  :: ycenter, ywindmin, ywindmax, edge, stressxTS, stressyTS, f
! phi, omega, f

  udif=0.;vdif=0.;
  f=2*OMEGA*sin(phi0deg*PI/180)
  !print *, f, dtf, step, dtime_dim, dtime, time_nondim, time_seconds
  !print *, phi, omega, f		
 
  !do i=1,NI
  !      print *,xc(i)
  !end do	 
 
  !do j=1,NJ
  !      print *,yc(j)
  !end do

  amplitude=0.05d0      
  
!*******************************************
! COMPUTATION OF THE WIND STRESS: NIHAR PAUL
!*******************************************
 
  ycenter=0.5*(yc(NJ)+yc(1));
  ywindmin=48.d0;
  ywindmax=yc(NJ)-48.d0;
  edge=0.06; 

  
  
  if (step .lt. 1000 .or. step .gt. 1600) then
                stressxTS=0.d0
                stressyTS=0.d0
  !             print *, stressxTS, stressyTS

  else
                stressxTS=amplitude*sin(f*dtime_dim*(step-1000))*0.5*(1-cos(2*PI*(step-1000)/600))
                stressyTS=amplitude*cos(f*dtime_dim*(step-1000))*0.5*(1-cos(2*PI*(step-1000)/600))
!		stressxTS=0.05d0
!		stressyTS=0.05d0
  !             print *, stressxTS, stressyTS
  endif

  do j=1,NJ
        if (yc(j) .lt.  yc(NJ/2)) then
                stressprofile(j)= 0.5*(tanh(edge*(yc(j)-ywindmin)*PI)+1.d0)
        else 
                stressprofile(j)=-0.5*(tanh(edge*(yc(j)-ywindmax)*PI)-1.d0)
        endif
        !print *, stressprofile(j)
  end do

  do j=1,NJ
        do i=1,NI
                stress_top_x(i,j)=stressxTS*stressprofile(j)
                stress_top_y(i,j)=stressyTS*stressprofile(j)
        end do
  end do
  !print *, stressprofile(48)
  !print *, stressxTS, stressyTS	
  !print *, stress_top_x(NI/2,NJ/2), stress_top_y(NI/2,NJ/2)
 
!****************************************
! COMPUTATION OF THE SOURCE TERM
!****************************************

  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 

  do j=1,NJ 
    do i=1,NI 
      stress_top(i,j) = sqrt(stress_top_x(i,j)*stress_top_x(i,j)+ stress_top_y(i,j)*stress_top_y(i,j))
      rhoinv = 1.d0/rho(i,j,NK) 
      !rhoinv = 1.d0/R0 
      Kdudzt= stress_top_x(i,j)*rhoinv*fact 
      Kdvdzt= stress_top_y(i,j)*rhoinv*fact 

      udif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdudzt   
      vdif(i,j,NK)= fac*Jac(i,j,NK)*wz(i,j,NK)*Kdvdzt   
    
    end do ! i
  end do ! j
                                                                        
return 
end subroutine wind_stress                                                                        
                                                                
