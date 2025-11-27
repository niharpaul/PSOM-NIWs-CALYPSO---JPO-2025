NFDIR='/usr/local'
gfortran ini_st.f90 -o exec_ini_st -I${NFDIR}/include -L${NFDIR}/lib -lnetcdff
