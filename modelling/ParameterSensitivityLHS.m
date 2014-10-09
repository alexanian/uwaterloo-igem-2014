function [ SensitivityCoeff, Errors ] = ...
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
        SampleCount = 1000;
    end
    
    PCount = size(ParameterBounds, 1);
    
    % Generate uniformly distributed parameters in [0,1] that
    % are independent of one another. 
    Samples = lhsdesign( SampleCount, PCount );
    for i = 1:PCount
        % Scale parameters using their ranges
        Samples(:,i) = ( ParameterBounds(i,2) - ParameterBounds(i, 1) ) * ...
            ( Samples(:,i) ) + ParameterBounds(i,1);
    end

    
    % Simulate DE for all parameter sets
    SampledErrors = zeros([SampleCount 1]);
    OutputFunc = @(Y) OutputSystem(Y);
    for i = 1:SampleCount
        DE = CreateParameterlessDE(DESystem, Samples(i,:));
        [~, Ys] = ode45( DE, TMeas, Y0 );
        
        % Calculate Error between Simulation and Measured Result
        SampledErrors(i) = LeastSquaresError( YMeas, OutputFunc(Ys) );
    end

    Errors = SampledErrors;

    % Calculate Threshold for deciding between
    % Acceptable and Unacceptable.
    Threshold = mean(SampledErrors);

    min(Errors)
    max(Errors)
    
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
        
        % Get ECDFs
        [f1, x1] = ecdf(Acceptable);
        [f2, x2] = ecdf(Unacceptable);
        
        % Plot difference of ECDFs (to observe sensitivity values across
        % the range of parameters)
        clf;
        range = linspace(ParameterBounds(i,1), ParameterBounds(i,2), SampleCount);
        semilogx(range, ecdfdiff(f1, x1, f2, x2, range));
        axis([-inf,inf,0,1]);
        
        % Save Difference in ECDF
        I = getframe(gcf);
        imwrite(I.cdata, ...
            sprintf('Plots/Sensitivity/System_Sensitivity_Difference_Parameter_%d.bmp',i));
        close(gcf)
        
    end
    
    figure('Position', [0 0 1024 1024]);
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
        
        % Plot ECDF of acceptable/unacceptable Parameter values
        subplot(4,4,i);
        ecdf(Acceptable);
        hold on;
        ecdf(Unacceptable);
        h = get(gca,'children');
        set(h(2), 'LineStyle', ':');
        title(sprintf('Parameter %d', i));
        hold off;
        
        % Return Sensitivity coefficient
        [~,~,SensitivityCoeff(i)] = kstest2(Acceptable, Unacceptable);
    end
    
    subplot(4,4,PCount+1);
    axis off;
    legend('Acceptable', 'Unacceptable');
    
    % Save CDF
    I = getframe(gcf);
    imwrite(I.cdata, ...
        sprintf('Plots/Sensitivity/System_Sensitivities.bmp'));
    close(gcf)
end

% Calculates the difference between two ECDFs
function FDiff = ecdfdiff( ecdf1, x1, ecdf2, x2, fullrange )
    c1 = 1;
    c2 = 1;
    c0 = 1;
    v1 = 0;
    v2 = 0;
    FDiff = zeros(size(fullrange));
    for x = fullrange
        while( x1(c1) < x )
            v1 = ecdf1(c1);
            c1 = c1 + 1;
            if( c1 > length(x1) )
                v1 = 1;
                c1 = length(x1);
                break;
            end
        end
        while( x2(c2) < x )
            v2 = ecdf2(c2);
            c2 = c2 + 1;
            if( c2 > length(x2) )
                v2 = 1;
                c2 = length(x2);
                break;
            end
        end
        
        FDiff(c0) = abs(v1 - v2);
        c0 = c0 + 1;
    end
end

% Calculates R Squared Error
function E = LeastSquaresError( YMeas, YObs )  
    %YM = mean(YObs);
    YRes = YObs - YMeas;
    %YTot = YObs - YM * ones(size(YObs));
    
    YRes = YRes .* YRes;
    %YTot = YTot .* YTot;
    
    E = sqrt(sum(YRes));
end

% Helper Function to Wrap a Parameterized DE System
function DE = CreateParameterlessDE( DEP, P )
    DE = @(T,Y) DEP(T,Y,P);
end