function dS = CRISPRi_System( t, s, p )
    % Production and degradation terms
    alpha_mrnaC = p(1);
    alpha_mrnaY = p(2);
    alpha_Rg = p(3);
    alpha0_mrnaY = 0;%p(4);
    beta_C = p(5);
    gamma_C = p(6);
    gamma_mrnaC = p(7);
    gamma_Rg = p(8);
    gamma_mrnaY = 0.1734;%p(9);

    % Reaction rate constants
    k_minus = p(10);
    k_plus = p(11);
    KA = p(13);
    n = p(14);
    
    % Differential terms
    mrnaC=s(1); % dCas9 mRNA
    C=s(2); % dCas9
    Rg=s(3); % sgRNA
    mrnaY=s(4); % mRNA to produce YFP

    % Non-differential terms from QSSA
    K = k_plus/k_minus;
    K_a = KA;
    CRg = K*C*Rg;

    mrnaC_dt=alpha_mrnaC - gamma_mrnaC*mrnaC;
    C_dt= beta_C*mrnaC - gamma_C*C + k_minus*CRg - k_plus*C*Rg;
    Rg_dt=alpha_Rg - gamma_Rg*Rg + k_minus*CRg - k_plus*C*Rg;
    mrnaY_dt=alpha_mrnaY/(1 + (CRg/K_a)^n) - gamma_mrnaY*mrnaY + alpha0_mrnaY;

    dS=[mrnaC_dt, C_dt, Rg_dt, mrnaY_dt]';
end

