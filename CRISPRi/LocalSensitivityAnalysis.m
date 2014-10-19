function ParameterRelSensitivity = ...
    LocalSensitivityAnalysis( DESystem, OutputSystem, Parameters, Y0, Tol )
    
    if nargin < 5
        Tol = 0.05;
    end
    
    DE = @(T,S) DESystem(T, S, Parameters);
    T = [0 ; 1000];
    [~, Y1] = ode45( DE, T, Y0 );
    Y1 = OutputSystem(Y1);
    
    ParameterRelSensitivity = zeros(size(Parameters));
    for p = 1:length(Parameters)
        Z = ones(size(Parameters));
        Z(p) = Tol + 1;
        DEP5 = @(T,S) DESystem(T, S, Parameters .* Z);
        
        DP = Parameters(p) * Tol;

        [~, Y2] = ode45( DEP5, T, Y0 );

        Y2 = OutputSystem(Y2);
        ParameterRelSensitivity(p) = ( Parameters(p) * ( Y2(length(Y2)) - ...
            Y1(length(Y1)) ) ) / ( DP * Y1(length(Y1)) );
    end
end