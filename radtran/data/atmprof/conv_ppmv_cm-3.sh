awk '{printf("%E\n", $11*1e-6*$2*1e3/($3*8.31*1e7)*6.02e23)}' tmp > rho_N2O

