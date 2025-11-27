subroutine ini_st 

	! --------------------                                              
	 USE header
	include 'netcdf.inc'

  	! This is the name of the data file we will read.
        character (len = *), parameter :: FILE_NAME = "/Users/niharpaul/psom/code/wiggle/src/psom_st_ic_1km.nc"

        ! We are reading 3D data, a 6 x 12 grid.
        integer, parameter :: NX=50, NY=98, NZ=34
        real(kind=8) :: s1(NX,NY,NZ), T1(NX,NY,NZ)

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid1, varid2

        ! Loop indexes, and error handling.
        integer :: i, j

        ! Create some pretend data. If this wasn't an example program, we
        ! would have some real data to write, for example, model output.
        ! do x = 1, NX
        !   do y = 1, NY
        !      data_out(y, x) = (x - 1) * NY + (y - 1)
        !   end do
        ! end do

        ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
        ! the file.
        call check(nf_open(FILE_NAME, NF_NOWRITE, ncid))

        ! Get the varid of the data variable, based on its name.
        call check(nf_inq_varid(ncid, "s", varid1))
        call check(nf_inq_varid(ncid, "T", varid2))

        ! Read the data.
        print*, "Reading Netcdf"
        print*, FILE_NAME

	call check(nf_get_var_double(ncid, varid1, s1))

        call check(nf_get_var_double(ncid, varid2, T1))
	
	call check(nf_close(ncid))

        print*,"*** SUCCESS reading example file: ", FILE_NAME, "! "
        print*, size(s1,DIM=1)
	print*, size(s1,DIM=2)
	print*, size(s1,DIM=3)
	print*, kind(s1)

	do j = 1, 1
           do i = 1, NX
              print*,s1(i,j,1)
           end do
        end do

	do j = 1, 1
           do i = 1, NX
              print*,T1(i,j,1)
           end do
        end do
 
	s(:,:,:,0)=s1
	T(:,:,:,0)=T1
        

	print*, size(s,DIM=1)
	print*, size(s,DIM=2)
	print*, size(s,DIM=3)  
	print*, "The value of salinity is:"
        print*,s(0,0,0,0)
        
	print*, "The value of temperature is:"
	print*,T(0,0,0,0)

	print*, kind(s)
	return

contains 
subroutine check(status)
		integer, intent (in) :: status
		integer :: trim
                if(status /= nf_noerr) then
                        print *, trim(nf_strerror(status))
                        stop "Stopped"
                end if
end subroutine check

end subroutine ini_st                                           
