function Error = Simulate_CRISPRi_Error( t_min, YD , Y0, Parameters )

    DE = @(T,S) CRISPRi_System( T, S, Parameters );
    [~, YDExperiment] = ode23( DE, t_min, Y0 );
    YDOut = CRISPRi_System_Output(YDExperiment);
    
    E = abs(YDOut - YD);
    E = E .* E;
    Error = sum(E);
end