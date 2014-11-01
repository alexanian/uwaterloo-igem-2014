function CRISPRi_Model
% CRISPRi_Model Wrapper for CRISPRi modelling/graphing & parameter estimation
%
% The model is based on two main reactions:
%
%   C + Rg <-k_minus, k_plus-> B
%   B <-k_2, k_1-> B_bound
%
% Where C is dCas9 protein, Rg is sgRNA and B is the dCas9-sgRnA complex
%
% The wrapper calls many of the scripts within this directory, including:
%
%   CRISPRi_DE_System.m:            model function to be used with ode45
%   CRISPRi_Plot.m:                 plots graphs in a consistent way
%   CRISPRi_Simulate_Time_Series.m: time series data based on Qi et al.
%   CRISPRi_Simulate_Error:         computes error (model vs. time series)
%   
%   See also CRISPRi_Literature_Parameters.m for initial parameters used in
%   the model & output data 35_fold_time_series.mat,6_fold_time_series.mat

    % Assign parameter values
    clear all
    
    % Load model parameters from the literature:
    load('CRISPRi_Literature_Parameters.mat');
    parameters_0 =      [alpha_mrnaC alpha_mrnaY alpha_Rg ...
                        beta_C gamma_C gamma_B gamma_mrnaC ...
                        gamma_Rg gamma_mrnaY k_minus k_plus K_a ...
                        n mrnaY_0 foldRepression];
    % Explanation of parameters:
    % mRNA production terms: alpha_mrnaC, alpha_mrnaY, alpha_Rg
    % mRNA degradation terms: gamma_mrnaC, gamma_mrnaY, gamma_Rg
    % dCas9 protein production: beta_C
    % protein degradation: gamma_C, gamma_B
    % reaction rates: k_minus, k_plus, K_a (=k_1/k_2)
    % hill constant: n
    % initial YFP mRNA: mrnaY_0
    % final n-fold repression: foldRepression

    % Set ODE simulation parameters
    CRISPRiODE=@CRISPRi_DE_System;
    options=odeset('Refine', 6);
    Tend=300;

    % Set initial condition: state = [mrnaC, C, Rg, mrnaY, mrnaY_noRepression]
    % Assume sgRNA starts at the same concentration as mrnaY as they are
    % under the same promoter
    S0=[0, 0, mrnaY_0, mrnaY_0, mrnaY_0]';

    % Run simulation
    [t,S]=ode45(CRISPRiODE, [0,Tend], S0, options, parameters_0);

    % Algebraic term from QSSA
    B = k_plus*S(:,2).*S(:,3)/(k_minus+gamma_B);

    % Plot initial results
    CRISPRi_Plot(t, S(:,1), S(:,2), S(:,3), B, S(:,4))
    title('Initial System Model of CRISPRi Repression')
    ylim([0 2])
    CRISPRi_Plot(t, S(:,1), S(:,2), S(:,3), B, S(:,5))
    title('Initial System Model without CRISPRi-YFP Interaction')
    ylim([0 2])

    % fmincon for improved mRNA production rates
    t_min = 1:300;
    
    % Expected YD (time series data to compare model to)
    YD_6 = CRISPRi_Simulate_Time_Series( t_min, mrnaY_0, 6);
    YD_35 = CRISPRi_Simulate_Time_Series( t_min, mrnaY_0, 35);
    figure; hold on
    plot(t_min, YD_6, '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    plot(t_min, YD_35, '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    hold off;

    % Initial parameter values and bounds.
    % note: only mRNA production rates are varied by fmincon
    parameters_0 =      [alpha_mrnaC alpha_mrnaY alpha_Rg ...
                        beta_C gamma_C gamma_B gamma_mrnaC ...
                        gamma_Rg gamma_mrnaY k_minus k_plus K_a ...
                        n mrnaY_0 foldRepression];
    parameters_upper =  [alpha_mrnaC*10^-3 alpha_mrnaY*10^-3 alpha_Rg*10^-3 ...
                        beta_C gamma_C gamma_B gamma_mrnaC ...
                        gamma_Rg gamma_mrnaY k_minus k_plus K_a ...
                        n mrnaY_0 foldRepression];
    parameters_lower =  [alpha_mrnaC*10^3 alpha_mrnaY*10^3 alpha_Rg*10^3 ...
                        beta_C gamma_C gamma_B gamma_mrnaC ...
                        gamma_Rg gamma_mrnaY k_minus k_plus K_a ...
                        n mrnaY_0 foldRepression];
                    
    presult_35 = fmincon(@(p) CRISPRi_Simulate_Error(t_min, YD_35, S0, p),...
        parameters_0, [], [], [], [], parameters_upper, parameters_lower );
    
    % Update the foldRepression to look at the worst-case repression
    foldRepression = 6;
    parameters_0(15) = foldRepression;
    parameters_upper(15) = foldRepression;
    parameters_lower(15) = foldRepression;
    
    presult_6 = fmincon(@(p) CRISPRi_Simulate_Error(t_min, YD_6, S0, p),...
        parameters_0, [], [], [], [], parameters_upper, parameters_lower )

    % Plot results for 35-fold repression
    [t_35,S_35]=ode45(CRISPRiODE, [0,Tend], S0, options, presult_35);
    S0=[0, 0, S_35(end,5), S_35(end,5), S_35(end,5)]';

    % Running a second time beginning at steady-state
    [t_35,S_35]=ode45(CRISPRiODE, [0,Tend], S0, options, presult_35);
    B = k_plus*S_35(:,2).*S_35(:,3)/(k_minus+gamma_B);

    CRISPRi_Plot(t_35, S_35(:,1), S_35(:,2), S_35(:,3), B, S_35(:,4))
    title('Adjusted System Model of CRISPRi Repression (35-fold)')
    ylim([0 2.1])
    plot(t_35, S_35(:,5), 'color', [0.98 0.93 0.36], 'Linewidth', 1.5)
    legend('dCas9 mRNA', 'dCas9', 'sgRNA', 'dCas9-sgRNA Complex',...
    'YFP mRNA', 'unrepressed YFP mRNA', 'Location','East')     
    save('Data_35FoldTimeSeries.mat', 't_35', 'S_35', 'B')

    
    % Plot results for 6-fold repression - start sgRNA and YFP at same
    % steady-state
    [t_6,S_6]=ode45(CRISPRiODE, [0,Tend], S0, options, presult_6);
    S0=[0, 0, S_6(end,5), S_6(end,5), S_6(end,5)]';

    % Running a second time beginning at steady-state
    [t_6,S_6]=ode45(CRISPRiODE, [0,Tend], S0, options, presult_6);
    B = k_plus*S_6(:,2).*S_6(:,3)/(k_minus+gamma_B);

    CRISPRi_Plot(t_6, S_6(:,1), S_6(:,2), S_6(:,3), B, S_6(:,4))
    title('Adjusted System Model of CRISPRi Repression (6-fold)')
    plot(t_6, S_6(:,5), 'color', [0.98 0.93 0.36], 'Linewidth', 1.5)
    legend('dCas9 mRNA', 'dCas9', 'sgRNA', 'dCas9-sgRNA Complex',...
    'YFP mRNA', 'unrepressed YFP mRNA', 'Location','East')
    ylim([0 2.1])    
    save('Data_6FoldTimeSeries.mat', 't_6', 'S_6', 'B')
    
    S_6(end,5)

    % Final plot comparing simulated time-series to fmincon result
    figure; hold on
    plot(t_6, CRISPRi_Simulate_Time_Series( t_6, mrnaY_0, 6 ), ...
        '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    plot(t_35, CRISPRi_Simulate_Time_Series( t_35, mrnaY_0, 35 ), ...
        '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    plot(t_6, S_6(:,4), 'color', [0.98 0.93 0.36], 'Linewidth', 3)
    plot(t_35, S_35(:,4), 'color', [0.98 0.93 0.36], 'Linewidth', 3)
    legend('6-Fold Time Series', '35-Fold Time Series', ...
        '6-Fold Model', '35-Fold Model');
    hold off

end