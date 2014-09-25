function SensitivityCoeff = ...
    ParameterSensitivityLHS( DESystem, OutputSystem, TMeas, YMeas, ParameterBounds, Y0, SampleCount )
% ParameterSensitivityLHS  Estimate Parameter Sensitivity for a Differential
% Equation System using Latin Hypercube Sampling
%
%       DESystem        
%                       Function with signature DY = ( T, Y, P )
%                       where (T,Y) is a point in the solution and P is the 
%                       (unknown) parameter vector that needs to have
%                       sensitivity determined.
%
%       OutputSystem    
%                       Takes the output vector from a solver and returns
%                       the output value that is associated with the 
%                       measurement.
%
%       TMeas,YMeas     
%                       Vector of measured values YMeas at timepoints
%                       TMeas.
%
%       ParameterBounds 
%                       A vector with dimensions [P,2] where P is the
%                       number of parameters and the 2 values represent the
%                       lower and upper bounds of the parameters.
%
%       Y0
%                       Initial Conditions for the DE (probably best if
%                       this aligns well with your measured data!)
%
%       SampleCount
%                       Number of sample simulations. Default is 10,000
%   Returns:
%       SensitivityCoeff
%                       A vector associating a parameter (index) to 
%                       a given value representing sensitivity. Higher
%                       values mean higher sensitivity.
%                       
%   See also lhsdesign, ode23.

    % Default Sample COunt
    if nargin < 7
        SampleCount = 10000;
    end
    
    PCount = size(ParameterBounds, 1);
    SampledErrors = zeros([SampleCount 1]);
    
    % Generate uniformly distributed parameters in [0,1] that
    % are independent of one another. 
    Samples = lhsdesign( SampleCount, PCount );
    for i = 1:PCount
        % Scale parameters using their ranges
        Samples(:,i) = ParameterBounds(i,2) * ...
            ( Samples(:,i) + ParameterBounds(i,1) );
    end

    % Simulate DE for all parameter sets
    parfor( i = 1:SampleCount, 20 )
        DE = CreateParameterlessDE(DESystem, Samples(i,:));
        [~, Ys] = ode23( DE, TMeas, Y0 );

        % Calculate Error between Simulation and Measured Result
        Ys = YMeas - OutputSystem(Ys);
        Ys = Ys .* Ys;
        SampledErrors(i) = sum(Ys);
    end

    % Calculate Threshold for deciding between
    % Acceptable and Unacceptable.
    Threshold = mean(SampledErrors);

    SensitivityCoeff = zeros([PCount 1]);
    % Go through each parameter and make a plot of their CDFS
    for i = 1:PCount

        % Partition different values of the parameter
        % into whether or not its simulation is acceptable or not
        Acceptable = zeros([SampleCount 1]);
        Unacceptable = zeros([SampleCount 1]);
        AC = 1;
        UC = 1;
        for s = 1:SampleCount
            if( SampledErrors(s) < Threshold )
                Acceptable(AC) = Samples(s,i);
                AC = AC + 1;
            else
                Unacceptable(UC) = Samples(s,i);
                UC = UC + 1;
            end
        end

        % Trim off excess parts of the vector...really doing this
        % for optimization purposes...
        Acceptable = Acceptable(1:(AC-1));
        Unacceptable = Unacceptable(1:(UC-1));
        
        % Plot CDF of Parameter value
        clf
        ecdf(Acceptable)
        hold on
        ecdf(Unacceptable)
        
        % Save CDF
        I = getframe(gcf);
        imwrite(I.cdata, ...
            sprintf('System_Sensitivity_Parameter_%d.bmp',i));
        close(gcf)
        
        % Get Sensitivity
        [~,~,SensitivityCoeff(i)] = kstest2(Acceptable, Unacceptable);
    end
end

% Helper Function to Wrap a Parameterized DE System
function DE = CreateParameterlessDE( DEP, P )
    DE = @(T,Y) DEP(T,Y,P);
end