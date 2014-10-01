function dS = CRISPRi_System( t, s, p )
    % Production and degradation terms
    alpha_mrnaC = p(1);
    alpha_mrnaY = p(2);
    alpha_Rg = p(3);
    alpha0_mrnaY = p(4);
    beta_C = p(5);
    gamma_C = p(6);
    gamma_mrnaC = p(7);
    gamma_Rg = p(8);
    gamma_mrnaY = p(9);

    % Reaction rate constants
    K = p(10);
    K_a = p(11);
    n = p(12);
    
    % Differential terms
    mrnaC=s(1); % dCas9 mRNA
    C=s(2); % dCas9
    Rg=s(3); % sgRNA
    mrnaY=s(4); % mRNA to produce YFP

    % Non-differential terms from QSSA
    CRg = K*C*Rg;

    % Actual System
    mrnaC_dt=alpha_mrnaC - gamma_mrnaC*mrnaC;
    C_dt= beta_C*mrnaC - gamma_C*C;
    Rg_dt=alpha_Rg - gamma_Rg*Rg;
    mrnaY_dt=alpha_mrnaY/(1 + (CRg/K_a)^n) - gamma_mrnaY*mrnaY + alpha0_mrnaY;

    dS=[mrnaC_dt, C_dt, Rg_dt, mrnaY_dt]';
end

