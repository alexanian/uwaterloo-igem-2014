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
    alpha_mrnaC=0.0011;%0.05;
    alpha_mrnaY=0.05;%0.05;
    alpha_Rg=0.001;%0.0011;
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
    mrnaY_0 = 0.29;
    foldRepression = 6;

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
    title('1');
    CRISPRi_Plot(t, S(:,1), S(:,2), S(:,3), B, S(:,5))
    title('2');

    % fmincon for improved mRNA production rates
    t_min = 1:270;
    Y0 = 5;

    % Expected YD (time series data to compare model to)
    YD_6 = Simulate_CRISPRi_Time_Series( t_min, Y0, 6);
    YD_35 = Simulate_CRISPRi_Time_Series( t_min, Y0, 35);

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

    S0 = [0, 0, 0, Y0, Y0]'; %state = [mrnaC, C, Rg, mrnaY, mrnaY_noRepression]
    presult = fmincon(@(p) Simulate_CRISPRi_Error(t_min, YD_35, S0, p),...
        parameters_0, [], [], [], [], parameters_upper, parameters_lower )
    Simulate_CRISPRi_Error(t_min, YD_35, S0, presult)
    
    % Plotzzz
    [t,S]=ode45(CRISPRiODE, [0,Tend], S0, options, presult);
    B = k_plus*S(:,2).*S(:,3)/(k_minus+gamma_B); % Algebraic term from QSSA

    CRISPRi_Plot(t, S(:,1), S(:,2), S(:,3), B, S(:,4))
    title('3');
    CRISPRi_Plot(t, S(:,1), S(:,2), S(:,3), B, S(:,5))
    title('4');

end