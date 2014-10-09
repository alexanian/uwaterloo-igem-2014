function Error = Simulate_CRISPRi_Error( TD, YD , Y0, Parameters )
    
    DE = @(T,S) CRISPRi_System( T, S, Parameters );
    [~, YDExperiment] = ode45( DE, TD, Y0 );
    YDOut = CRISPRi_system_output(YDExperiment);
    
    E = abs(YDOut - YD);
    E = E .* E;
    Error = sum(E);
end