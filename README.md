# PSOM-NIWs-CALYPSO---JPO-2025
The model code used to study the transport of passive tracers driven by Near-inertial Waves in a density front.
CDF refers to the CALYPSO Density Front. The initial front is generated in the file ini_st.f90, while the tracer is initialized using the profile specified in ini_tracer.f90. The wind-stress field can be prescribed through windstress.f90. A suitable Hanning window is applied to the NIW forcing to prevent spurious jumps in the simulation. The modified model source codes are located in cdf_00/src. 


