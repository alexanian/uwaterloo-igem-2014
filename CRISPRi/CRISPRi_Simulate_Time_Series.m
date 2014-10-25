function YD = Simulate_CRISPRi_Time_Series( t_min, Y0, FR )
% Define exponential using 35-minute half life from Qi et al. 2013 
    lambda = log(2)/35;
    
    YD = Y0*exp(-lambda*(t_min-8));
    YD(YD > Y0) = Y0; 
    YD(YD < Y0/FR) = Y0/FR;
    YD = YD';
end