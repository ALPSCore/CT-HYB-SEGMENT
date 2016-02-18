npts_tau = 256
python gen_sample_hybridization.py --gamma 0.2 --npts ${npts_tau}
alps_cthyb --N_ORBITALS 2 --U 4 --MU 2 --N_TAU ${npts_tau} --BETA 1 --DELTA delta_tau.dat --N_MEAS 1000 --THERMALIZATION 100 --SWEEPS 200000 --TEXT_OUTPUT 1
