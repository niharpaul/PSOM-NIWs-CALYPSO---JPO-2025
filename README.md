# PSOM-NIWs-CALYPSO---JPO-2025
The model code used to study the transport of passive tracers driven by Near-inertial Waves in a density front.
CDF refers to the CALYPSO Density Front. The initial front is generated in ini_st.f90, and the tracer is initialized using the profile specified in ini_tracer.f90. The wind-stress field can be prescribed through windstress.f90. A Hanning window is applied to the NIW forcing to avoid spurious jumps in the simulation. The modified model source codes are located in cdf_00/src. Certain inputs for running the model must be hardcoded. The model timeframe, data formats, and general input parameters are specified in cdf_00/namelist_cdf_00.


