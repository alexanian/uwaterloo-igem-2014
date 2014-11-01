function ParameterRelSensitivity = ...
    LocalSensitivityAnalysis( DESystem, System_YFP_Output, Parameters, Y0, Tol )
    
    if nargin < 5
        Tol = 0.05;
    end
    
    DE = @(T,S) DESystem(T, S, Parameters);
    T = [0 ; 1000];
    [~, Y1] = ode45( DE, T, Y0 );
    YFP1 = System_YFP_Output(Y1);
    
    ParameterRelSensitivity = zeros(size(Parameters));
    for p = 1:length(Parameters)
        Z = ones(size(Parameters));
        Z(p) = Tol + 1;
        DEP5 = @(T,S) DESystem(T, S, Parameters .* Z);
        
        DP = Parameters(p) * Tol;

        [~, Y2] = ode45( DEP5, T, Y0 );

        YFP2 = System_YFP_Output(Y2);
        ParameterRelSensitivity(p) = ( Parameters(p) * ( YFP2 - YFP1 ) )...
        / ( DP * YFP1 );
    end
end