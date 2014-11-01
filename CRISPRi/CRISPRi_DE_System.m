function dS = CRISPRi_DE_System( ~, s, p )
% CRISPRi_System Model of CRISPRi system to be used with ode45/ode23
    % Production and degradation terms
    alpha_mrnaC = p(1);
    alpha_mrnaY = p(2);
    alpha_Rg = p(3);
    beta_C = p(4);
    gamma_C = p(5);
    gamma_B = p(6);
    gamma_mrnaC = p(7);
    gamma_Rg = p(8);
    gamma_mrnaY = p(9);

    % Reaction rate constants
    k_minus = p(10);
    k_plus = p(11);
    K_a = p(12);
    n = p(13);
    
    % Maximal concentration parameters
    mrnaY_0 = p(14);
    foldRepression = p(15);

    % Differential terms
    mrnaC=s(1); % dCas9 mRNA
    C=s(2); % dCas9
    Rg=s(3); % sgRNA
    mrnaY=s(4); % mRNA to produce YFP
    mrnaY_noR=s(5); % mRNA dynamics without repression

    % Non-differential terms from QSSA
    B = k_plus*C*Rg/(k_minus+gamma_B);

    % Differential System
    mrnaC_dt =      alpha_mrnaC - gamma_mrnaC*mrnaC;
    C_dt =          beta_C*mrnaC - gamma_C*C + k_minus*B - k_plus*C*Rg;
    Rg_dt =         alpha_Rg - gamma_Rg*Rg + k_minus*B - k_plus*C*Rg;
    
    alpha_ratio =   gamma_mrnaY*mrnaY_0/(alpha_mrnaY*foldRepression);
    
    mrnaY_dt =      alpha_mrnaY*(alpha_ratio + ...
                    (1-alpha_ratio)/(1 + (B/K_a)^n)) - gamma_mrnaY*mrnaY;
    mrnaY_noR_dt =  alpha_mrnaY - gamma_mrnaY*mrnaY_noR;
    
    dS=[mrnaC_dt, C_dt, Rg_dt, mrnaY_dt, mrnaY_noR_dt]';
end

