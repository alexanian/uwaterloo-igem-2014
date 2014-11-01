function Error = CRISPRi_Simulate_Error( t_min, YD , S0, Parameters )
% Find squared error comparing simulated time series with model output
    DE = @(T,S) CRISPRi_DE_System( T, S, Parameters );
    [~, YDExperiment] = ode23( DE, t_min, S0 );
    YDOut = CRISPRi_System_YFP_Output(YDExperiment);
    
    E = abs(YDOut - YD);
    E = E .* E;
    Error = sum(E);
end