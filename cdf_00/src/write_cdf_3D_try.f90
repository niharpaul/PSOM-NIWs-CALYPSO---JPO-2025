subroutine write_cdf_3D(stepl,n)

!----------------------------------------------------------             
! 3D output routine.
! It writes both center and face values to be written in .cdf files.
!----------------------------------------------------------             

    USE header, only : NI,NJ,NK,xc,yc,zc,p,h,consump,T,s,rho,Tr,u,v,w,vor,pv,uf,vf,wf,Kz,conv,con100,nconsume,dirout,rc_kind,ntr
#include "netcdf.inc"                                                                        

    integer :: n,stepl,i,j,k,it

    character (len = 100) :: FILE_NAME

    integer, parameter :: NJdir = NJ+2, NIdir = NI+2, NKdir = NK+2, ntrdir = ntr
    ! Dimensions
    integer :: xdir_dimid, ydir_dimid, zdir_dimid, ntr_dimid
    integer :: xdir_varid, ydir_varid, zdir_varid
    ! 2D variable
    integer :: h_varid
    ! 3D variable
    integer :: p_varid, T_varid, S_varid, rho_varid, u_varid, v_varid, w_varid, vor_varid, pv_varid, c_varid, c100_varid, uf_varid, vf_varid, wf_varid, Kz_varid
    ! 4D variable
    integer :: Tr_varid, cp_varid

    integer :: dimids2(2)
    integer :: dimids3(3)
    integer :: dimids4(4)
    integer :: ncid

    REAL(kind=rc_kind) :: Trwrite(0:NI+1,0:NJ+1,0:NK+1,ntr)
    REAL(kind=rc_kind) :: z(0:NI+1,0:NJ+1,0:NK+1) 
    REAL(kind=rc_kind) :: zf(0:NI+1,0:NJ+1,0:NK+1) 
    REAL(kind=rc_kind) :: xf(0:NI+1,0:NJ+1,0:NK+1) 
    REAL(kind=rc_kind) :: yf(0:NI+1,0:NJ+1,0:NK+1) 
    REAL(kind=rc_kind) ::  rcode
    
    ! For the sake of better plots 
    do k=0,NK+1 
        do j=0,NJ+1 
            do i=0,NI+1 
                z(i,j,k)= zc(i,j,k) 
                do it=1,ntr 
                    Trwrite(i,j,k,it)= Tr(it,i,j,k,n) 
                end do
            end do
        end do
    end do

    do j=0,NJ+1 
        do i=0,NI+1 
            z(i,j,NK+1)= 0.d0 ;   z(i,j,0)= 0.5*(z(i,j,0)+z(i,j,1)); 
            s(i,j,NK+1,n)= s(i,j,NK,n); s(i,j,0,n)= s(i,j,1,n); 
            T(i,j,NK+1,n)= T(i,j,NK,n); T(i,j,0,n)= T(i,j,1,n); 
            rho(i,j,NK+1)= rho(i,j,NK); rho(i,j,0)= rho(i,j,1); 
            u(i,j,NK+1,n)= u(i,j,NK,n); u(i,j,0,n)= u(i,j,1,n); 
            v(i,j,NK+1,n)= v(i,j,NK,n); v(i,j,0,n)= v(i,j,1,n); 
            w(i,j,NK+1,n)= w(i,j,NK,n); w(i,j,0,n)= w(i,j,1,n); 
            vor(i,j,NK+1)= vor(i,j,NK); vor(i,j,0)= vor(i,j,1); 
            pv(i,j,NK+1) = pv(i,j,NK) ; pv(i,j,0) = pv(i,j,1) ; 
            do it=1,ntr 
                Trwrite(i,j,NK+1,it)= Trwrite(i,j,NK,it); Trwrite(i,j,0,it)= Trwrite(i,j,1,it); 
            end do
        end do
    end do
        
!    print*, zc

    ! CHECK UNITS OF EACH VARIABLES WRITTEN TO MAKE SURE IT MATCHES SI UNITS
    ! ==> scale xc yc and zc
    ! MAKE SURE CELL FACES VECTORS ARE IMPORTED AND WRITTEN
    ! ==> FIND XF and YF

    !---------------------------------------------------------
    !---------------------------------------------------------
    ! FULL FILE - CELL CENTER VARIABLES
    !---------------------------------------------------------
    !---------------------------------------------------------

    ! This is the name of the data file we will create.
    WRITE(FILE_NAME,'("full_",I5.5,".cdf")') stepl     ! Cell centers 

    ! Create the file. 
    ncid =  nccre(TRIM(dirout)//FILE_NAME,NCCLOB,rcode)

    !!!!!!!!!!!!!!!!!
    ! DIMENSIONS
    !!!!!!!!!!!!!!!!!

    ! Define the dimensions.
    xdir_dimid = ncddef(ncid, "i-dim", NIdir, rcode)
    ydir_dimid = ncddef(ncid, "j-dim", NJdir, rcode)
    zdir_dimid = ncddef(ncid, "k-dim", NKdir, rcode)
    ntr_dimid = ncddef(ncid, "ntr-dim", ntrdir, rcode)

    !!!!!!!!!!!!!!!!!
    ! COORDINATES
    !!!!!!!!!!!!!!!!!
    ! Define the coordinates. A varid is returned for each.
    xdir_varid = ncvdef(ncid, "xc", NCDOUBLE, 1, xdir_dimid, rcode)
    ydir_varid = ncvdef(ncid, "yc", NCDOUBLE, 1, ydir_dimid, rcode)
    zdir_varid = ncvdef(ncid, "zc", NCDOUBLE, 3, (/ xdir_dimid, ydir_dimid, zdir_dimid /), rcode)
    
    ! Assign attributes to coordinates.
    CALL NCAPTC (ncid, xdir_varid,'long_name', NCCHAR, 36,'cell center coordinate - x direction', rcode)
    CALL NCAPTC (ncid, xdir_varid,'units', NCCHAR, 5,'meter', rcode)
    CALL NCAPTC (ncid, ydir_varid,'long_name', NCCHAR, 36,'cell center coordinate - y direction', rcode)
    CALL NCAPTC (ncid, ydir_varid,'units', NCCHAR, 5,'meter', rcode)
    CALL NCAPTC (ncid, zdir_varid,'long_name', NCCHAR, 36,'cell center coordinate - z direction', rcode)
    CALL NCAPTC (ncid, zdir_varid,'units', NCCHAR, 5,'meter', rcode)        
    
    CALL ncendf(ncid,rcode)
    
    ! Write the coordinate data.
    CALL ncvpt(ncid,xdir_varid, 1, NIdir, xc, rcode)
    CALL ncvpt(ncid,ydir_varid, 1, NJdir, yc, rcode)
    CALL ncvpt(ncid,zdir_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), z*1000, rcode)

    !!!!!!!!!!!!!!!!!
    ! DATA
    !!!!!!!!!!!!!!!!!
    CALL NCREDF(ncid, rcode)
    ! Define variables. A varid is returned for each.
    ! 2D variables
    dimids2 = (/ xdir_dimid, ydir_dimid /)
    h_varid = ncvdef(ncid, "h", NCDOUBLE, 2, dimids2, rcode)

    ! 3D variables
    dimids3 = (/ xdir_dimid, ydir_dimid, zdir_dimid /)
    p_varid = ncvdef(ncid, "p", NCDOUBLE, 3, dimids3, rcode)
    T_varid = ncvdef(ncid, "T", NCDOUBLE, 3, dimids3, rcode)
    S_varid = ncvdef(ncid, "S", NCDOUBLE, 3, dimids3, rcode)
    rho_varid = ncvdef(ncid, "rho", NCDOUBLE, 3, dimids3, rcode)
    u_varid = ncvdef(ncid, "u", NCDOUBLE, 3, dimids3, rcode)
    v_varid = ncvdef(ncid, "v", NCDOUBLE, 3, dimids3, rcode)
    w_varid = ncvdef(ncid, "w", NCDOUBLE, 3, dimids3, rcode)
    vor_varid = ncvdef(ncid, "vor", NCDOUBLE, 3, dimids3, rcode)
    pv_varid = ncvdef(ncid, "pv", NCDOUBLE, 3, dimids3, rcode)
    !c_varid = ncvdef(ncid, "conv", NCDOUBLE, 3, dimids3, rcode)
    !c100_varid = ncvdef(ncid, "con100", NCDOUBLE, 3, dimids3, rcode)


    ! 4D variables
    dimids4 = (/ xdir_dimid, ydir_dimid, zdir_dimid, ntr_dimid /)
    Tr_varid = ncvdef(ncid, "Tr", NCDOUBLE, 4, dimids4, rcode)
    !dimids4 = (/ xdir_dimid, ydir_dimid, zdir_dimid, nconsume /)
    !cp_varid = ncvdef(ncid, "consump", NCDOUBLE, 4, dimids4, rcode)

    ! Assign attributes to variables.
    ! 2D variables
    CALL NCAPTC (ncid, h_varid,'standard_name', NCCHAR, 18,'sea_surface_height', rcode)
    CALL NCAPTC (ncid, h_varid,'unit', NCCHAR, 5,'meter', rcode)
    CALL NCAPTC (ncid, h_varid,'long_name', NCCHAR, 18,'Sea surface height', rcode)

    ! 3D variables
    CALL NCAPTC (ncid, p_varid,'standard_name', NCCHAR, 18,'sea_water_pressure', rcode)
    CALL NCAPTC (ncid, p_varid,'unit', NCCHAR, 4,'dbar', rcode)
    CALL NCAPTC (ncid, p_varid,'long_name', NCCHAR, 8,'Pressure', rcode)

    CALL NCAPTC (ncid, T_varid,'standard_name', NCCHAR, 21,'sea_water_temperature', rcode)
    CALL NCAPTC (ncid, T_varid,'unit', NCCHAR, 7,'Celsius', rcode)
    CALL NCAPTC (ncid, T_varid,'long_name', NCCHAR, 11,'Temperature', rcode)

    CALL NCAPTC (ncid, S_varid,'standard_name', NCCHAR, 18,'sea_water_salinity', rcode)
    CALL NCAPTC (ncid, S_varid,'unit', NCCHAR, 1,'', rcode)
    CALL NCAPTC (ncid, S_varid,'long_name', NCCHAR, 8,'Salinity', rcode)
    
    CALL NCAPTC (ncid, rho_varid,'standard_name', NCCHAR, 17,'sea_water_density', rcode)
    CALL NCAPTC (ncid, rho_varid,'unit', NCCHAR, 6,'kg m-3', rcode)
    CALL NCAPTC (ncid, rho_varid,'long_name', NCCHAR, 7,'Density', rcode)

    CALL NCAPTC (ncid, u_varid,'standard_name', NCCHAR, 27,'eastward_sea_water_velocity', rcode)
    CALL NCAPTC (ncid, u_varid,'unit', NCCHAR, 5,'m s-1', rcode)
    CALL NCAPTC (ncid, u_varid,'long_name', NCCHAR, 17,'Eastward velocity', rcode)

    CALL NCAPTC (ncid, v_varid,'standard_name', NCCHAR, 28,'northward_sea_water_velocity', rcode)
    CALL NCAPTC (ncid, v_varid,'unit', NCCHAR, 5,'m s-1', rcode)
    CALL NCAPTC (ncid, v_varid,'long_name', NCCHAR, 18,'Northward velocity', rcode)
    
    CALL NCAPTC (ncid, w_varid,'standard_name', NCCHAR, 25,'upward_sea_water_velocity', rcode)
    CALL NCAPTC (ncid, w_varid,'unit', NCCHAR, 5,'m s-1', rcode)
    CALL NCAPTC (ncid, w_varid,'long_name', NCCHAR, 17,'Vertical velocity', rcode)
    
    CALL NCAPTC (ncid, vor_varid,'standard_name', NCCHAR, 24,'ocean_relative_vorticity', rcode)
    CALL NCAPTC (ncid, vor_varid,'unit', NCCHAR, 3,'s-1', rcode)
    CALL NCAPTC (ncid, vor_varid,'long_name', NCCHAR, 18,'relative vorticity', rcode)
    
    CALL NCAPTC (ncid, pv_varid,'standard_name', NCCHAR, 34,'potential_vorticity_of_ocean_layer', rcode)
    CALL NCAPTC (ncid, pv_varid,'unit', NCCHAR, 7,'m-1 s-1', rcode)
    CALL NCAPTC (ncid, pv_varid,'long_name', NCCHAR, 19,'Potential vorticity', rcode)
    
    ! 4D variables
    CALL NCAPTC (ncid, Tr_varid,'unit', NCCHAR, 1,'', rcode)
    CALL NCAPTC (ncid, Tr_varid,'long_name', NCCHAR, 20,'Tracer concentration', rcode)
    
    CALL ncendf(ncid,rcode)
    
    ! Write the variable data.
    ! 2D variables
    CALL ncvpt(ncid,h_varid, (/ 1,1 /), (/ NIdir, NJdir /), h, rcode)

    ! 3D variables
    CALL ncvpt(ncid,p_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), p, rcode)
    CALL ncvpt(ncid,T_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), T, rcode)
    CALL ncvpt(ncid,S_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), S, rcode)
    CALL ncvpt(ncid,rho_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), rho, rcode)
    CALL ncvpt(ncid,u_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), u, rcode)
    CALL ncvpt(ncid,v_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), v, rcode)
    CALL ncvpt(ncid,w_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), w, rcode)
    CALL ncvpt(ncid,vor_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), vor, rcode)
    CALL ncvpt(ncid,pv_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), pv, rcode)
    CALL ncvpt(ncid,c_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), conv, rcode)
    CALL ncvpt(ncid,c100_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), con100, rcode)

    ! 4D variables
    CALL ncvpt(ncid,Tr_varid, (/ 1,1,1,1 /), (/ NIdir, NJdir, NKdir, ntr /), Trwrite, rcode)
    !CALL ncvpt(ncid,cp_varid, (/ 1,1,1,1 /), (/ NIdir, NJdir, NKdir, ntr /), consump, rcode)

    ! Close the file.
    CALL ncclos(ncid,rcode)
 
    ! If we got this far, everything worked as expected. Yipee! 
    print *,"*** SUCCESS writing full file output!"


    !---------------------------------------------------------
    !---------------------------------------------------------
    ! FACE FILE - CELL CENTER VARIABLES
    !---------------------------------------------------------
    !---------------------------------------------------------
    ! This is the name of the data file we will create.
    WRITE(FILE_NAME,'("face_",I5.5,".cdf")') stepl     ! Cell centers 

    ! Create the file. 
    ncid =  nccre(TRIM(dirout)//FILE_NAME,NCCLOB,rcode)

    !!!!!!!!!!!!!!!!!
    ! DIMENSIONS
    !!!!!!!!!!!!!!!!!

    ! Define the dimensions.
    xdir_dimid = ncddef(ncid, "i-dim", NIdir, rcode)
    ydir_dimid = ncddef(ncid, "j-dim", NJdir, rcode)
    zdir_dimid = ncddef(ncid, "k-dim", NKdir, rcode)
    
    !!!!!!!!!!!!!!!!!
    ! COORDINATES
    !!!!!!!!!!!!!!!!!
    ! Define the coordinates. A varid is returned for each.
    xdir_varid = ncvdef(ncid, "xf", NCDOUBLE, 1, xdir_dimid, rcode)
    ydir_varid = ncvdef(ncid, "yf", NCDOUBLE, 1, ydir_dimid, rcode)
    zdir_varid = ncvdef(ncid, "zf", NCDOUBLE, 3, (/ xdir_dimid, ydir_dimid, zdir_dimid /), rcode)

    ! Assign attributes to coordinates.
    CALL NCAPTC (ncid, xdir_varid,'long_name', NCCHAR, 34,'cell face coordinate - x direction', rcode)
    CALL NCAPTC (ncid, xdir_varid,'units', NCCHAR, 5,'meter', rcode)
    CALL NCAPTC (ncid, ydir_varid,'long_name', NCCHAR, 34,'cell face coordinate - y direction', rcode)
    CALL NCAPTC (ncid, ydir_varid,'units', NCCHAR, 5,'meter', rcode)
    CALL NCAPTC (ncid, zdir_varid,'long_name', NCCHAR, 34,'cell face coordinate - z direction', rcode)
    CALL NCAPTC (ncid, zdir_varid,'units', NCCHAR, 5,'meter', rcode)        
    
    CALL ncendf(ncid,rcode)
    
    ! Write the coordinate data.
!    CALL ncvpt(ncid,xdir_varid, 1, NIdir, xf, rcode)
!    CALL ncvpt(ncid,ydir_varid, 1, NJdir, yf, rcode)
!    CALL ncvpt(ncid,zdir_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), zf, rcode)

    !!!!!!!!!!!!!!!!!
    ! DATA
    !!!!!!!!!!!!!!!!!
    CALL NCREDF(ncid, rcode)
    
    ! Define variables. A varid is returned for each.
    ! 3D variables
    dimids3 = (/ xdir_dimid, ydir_dimid, zdir_dimid /)
    uf_varid = ncvdef(ncid, "uf", NCDOUBLE, 3, dimids3, rcode)
    vf_varid = ncvdef(ncid, "vf", NCDOUBLE, 3, dimids3, rcode)
    wf_varid = ncvdef(ncid, "wf", NCDOUBLE, 3, dimids3, rcode)
    Kz_varid = ncvdef(ncid, "Kz", NCDOUBLE, 3, dimids3, rcode)
    
    ! Assign attributes to variables.
    ! 3D variables
    CALL NCAPTC (ncid, uf_varid,'standard_name', NCCHAR, 1,'', rcode)
    CALL NCAPTC (ncid, uf_varid,'unit', NCCHAR, 9,'m s-1 m-2', rcode)
    CALL NCAPTC (ncid, uf_varid,'long_name', NCCHAR, 50,'Eastward velocity flux through the model cell face', rcode)
    CALL NCAPTC (ncid, uf_varid,'notes', NCCHAR, 34,'uf, vf, and wf are divergence-free', rcode)

    CALL NCAPTC (ncid, vf_varid,'standard_name', NCCHAR, 1,'', rcode)
    CALL NCAPTC (ncid, vf_varid,'unit', NCCHAR, 9,'m s-1 m-2', rcode)
    CALL NCAPTC (ncid, vf_varid,'long_name', NCCHAR, 51,'Northward velocity flux through the model cell face', rcode)
    CALL NCAPTC (ncid, vf_varid,'notes', NCCHAR, 34,'uf, vf, and wf are divergence-free', rcode)

    CALL NCAPTC (ncid, wf_varid,'standard_name', NCCHAR, 26,'ocean_vertical_diffusivity', rcode)
    CALL NCAPTC (ncid, wf_varid,'unit', NCCHAR, 6,'m2 s-1', rcode)
    CALL NCAPTC (ncid, wf_varid,'long_name', NCCHAR, 20,'Vertical Diffusivity', rcode)

    CALL ncendf(ncid,rcode)

    ! Write the variable data.
    ! 3D variables
    CALL ncvpt(ncid,uf_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), uf, rcode)
    CALL ncvpt(ncid,vf_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), vf, rcode)
    CALL ncvpt(ncid,wf_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), wf, rcode)
    CALL ncvpt(ncid,Kz_varid, (/ 1,1,1 /), (/ NIdir, NJdir, NKdir /), Kz, rcode)

    ! Close the file.
    CALL ncclos(ncid,rcode)
 
    ! If we got this far, everything worked as expected. Yipee! 
    print *,"*** SUCCESS writing face file output!"

    return
end