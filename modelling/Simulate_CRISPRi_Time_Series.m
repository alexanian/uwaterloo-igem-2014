function YD = Simulate_CRISPRi_Time_Series( t_min, Y0, FR )
    % Define exponential using data from Qi et al. 2013
    YD = Y0*exp(-(t_min-8)/35);
    YD(YD > Y0) = Y0; 
    YD(YD < Y0/FR) = Y0/FR;
    YD = YD';
end