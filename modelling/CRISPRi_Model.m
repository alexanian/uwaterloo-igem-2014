% Prepares graphs and declares parameters for the CRISPRi_System
% ODEs
% Two main reactions:
% C + Rg <-k_minus, k_plus-> CRg
% CRg + O <-k_2, k_1-> CRgO
% Where "O" I believe is a repressor for YFP

function CRISPRi_Model
    % Assign parameter values
    clear all
    
    % Production and degradation terms
    alpha_mrnaC=0.0011;
    alpha_mrnaY=0.0011;
    alpha_Rg=0.0011;
    beta_C=(0.0057 + 0.4797)/2; % midpoint of range
    gamma_C=5.6408*10^-4;
    gamma_B=5.6408*10^-4;
    gamma_mrnaC=0.1734;
    gamma_Rg=0.1734;
    gamma_mrnaY=0.1734;

    % Reaction rate constants
    k_minus=0.1; % midpoint of range
    k_plus=0.1; % midpoint of range
    K_a=0.28;
    n=2.5;

    % Maximal concentration parameters
    mrnaY_0 = 2;
    foldRepression = 35;

    % Set simulation parameters
    CRISPRiODE=@CRISPRi_System;
    options=odeset('Refine', 6);
    Tend=300;
    parameters_0 =      [alpha_mrnaC alpha_mrnaY alpha_Rg ...
                        beta_C gamma_C gamma_B gamma_mrnaC ...
                        gamma_Rg gamma_mrnaY k_minus k_plus K_a ...
                        n mrnaY_0 foldRepression];

    % Set initial condition: state = [mrnaC, C, Rg, mrnaY, mrnaY_noRepression]
    S0=[0, 0, 0, mrnaY_0, mrnaY_0]';

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
    Y0 = 2;

    % Expected YD (time series data to compare model to)
    YD_6 = Simulate_CRISPRi_Time_Series( t_min, Y0, 6);
    YD_35 = Simulate_CRISPRi_Time_Series( t_min, Y0, 35);
    figure; hold on
    plot(t_min, YD_6, '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    plot(t_min, YD_35, '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    hold off;

    % Initial parameter values and bounds; note that only mRNA production rates
    % are varied
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
    presult_35 = fmincon(@(p) Simulate_CRISPRi_Error(t_min, YD_35, S0, p),...
        parameters_0, [], [], [], [], parameters_upper, parameters_lower )
    
    foldRepression = 6;
    parameters_0(15) = foldRepression;
    parameters_upper(15) = foldRepression;
    parameters_lower(15) = foldRepression;
    presult_6 = fmincon(@(p) Simulate_CRISPRi_Error(t_min, YD_6, S0, p),...
        parameters_0, [], [], [], [], parameters_upper, parameters_lower )

    % Plot results for 6-fold repression
    [t_6,S_6]=ode45(CRISPRiODE, [0,Tend], S0, options, presult_6);
    B = k_plus*S_6(:,2).*S_6(:,3)/(k_minus+gamma_B);

    CRISPRi_Plot(t_6, S_6(:,1), S_6(:,2), S_6(:,3), B, S_6(:,4))
    title('Adjusted System Model of CRISPRi Repression (6-fold)')
    ylim([0 2])
    CRISPRi_Plot(t_6, S_6(:,1), S_6(:,2), S_6(:,3), B, S_6(:,5))
    title('Adjusted System Model without CRISPRi-YFP Interaction (6-fold)')
    ylim([0 2])

    % Plot results for 35-fold repression
    [t_35,S_35]=ode45(CRISPRiODE, [0,Tend], S0, options, presult_35);
    B = k_plus*S_35(:,2).*S_35(:,3)/(k_minus+gamma_B);

    CRISPRi_Plot(t_35, S_35(:,1), S_35(:,2), S_35(:,3), B, S_35(:,4))
    title('Adjusted System Model of CRISPRi Repression (35-fold)')
    ylim([0 2])
    CRISPRi_Plot(t_35, S_35(:,1), S_35(:,2), S_35(:,3), B, S_35(:,5))
    title('Adjusted System Model without CRISPRi-YFP Interaction (35-fold)')
    ylim([0 2])
    
    % Final plot comparing simulated time-series to fmincon result
    size(t_6)
    size(S_35)
    size(S_35)
    figure; hold on
    plot(t_6, Simulate_CRISPRi_Time_Series( t_6, Y0, 6 ), ...
        '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    plot(t_35, Simulate_CRISPRi_Time_Series( t_35, Y0, 35 ), ...
        '--', 'color', [0.93 0.73 0.36], 'Linewidth', 1.5)
    plot(t_6, S_6(:,4), 'color', [0.98 0.93 0.36], 'Linewidth', 3)
    plot(t_35, S_35(:,4), 'color', [0.98 0.93 0.36], 'Linewidth', 3)
    legend('6-Fold Time Series', '35-Fold Time Series', ...
        '6-Fold Model', '35-Fold Model');
    hold off

end