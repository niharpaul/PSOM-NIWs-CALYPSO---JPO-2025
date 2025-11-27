subroutine tracerinit(stepl) 
        USE header, only : Tr, zc, NI, NJ, NK, DL
        integer :: it,stepl, i, j, k
!       initializes tracer fields                                         
!       TRACER 1                                                          
!       =========                                                         
        it = 1 
        Tr = 0d0
        zc=zc*DL
        if (stepl.eq.1600) then
        !zm=zc(1:NI,1:NJ,2:NK+1)+zc(1:NI,1:NJ,0:K)
        do i=1,NI
                do j=1,NJ
                        do k=1,NK
                                Tr(it,i,j,k,0) =(zc(i,j,k)/(zc(i,j,NK)-zc(i,j,1)))-(zc(i,j,1)/(zc(i,j,NK)-zc(i,j,1)))
                                !print *, Tr(it,i,j,2,0)
                        end do
                        !Tr(it,i,j,0,0)=Tr(it,i,j,1,0)
                        !Tr(it,i,j,NK+1,0)=Tr(it,i,j,NK,0)
                end do
        end do
        !Tr(it,0,1:NJ,0:NK+1,0)=Tr(it,1,1:NJ,0:NK+1,0)
        !Tr(it,NI+1,1:NJ,0:NK+1,0)=Tr(it,NI,1:NJ,0:NK+1,0)
        !Tr(it,1:NI,0,0:NK+1,0)=Tr(it,1:NI,1,0:NK+1,0)
        !Tr(it,1:NI,NJ+1,0:NK+1,0)=Tr(it,1:NI,NJ,0:NK+1,0)
        endif

        !do k=1,NK
        !print *, zc(NI/2,NJ/2,k), Tr(it,NI/2,NJ/2,k,0)
        !end do
end subroutine tracerinit
