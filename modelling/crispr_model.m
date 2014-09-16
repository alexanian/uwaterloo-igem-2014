%Basic model for the CRISPR system
% Two main reactions:
% C + Rg <-k_minus, k_plus-> CRg
% CRg + O <-k_2, k_1-> CRgO
% Where "O" I believe is a repressor for YFP

function crispr_model

clear all

%declare model parameters
% Production and degradation terms
global alpha_mrnaC
global alpha_mrnaY
global alpha_Rg
global alpha0_mrnaY
global beta_C
global gamma_C
global gamma_mrnaC
global gamma_Rg
global gamma_mrnaY

% Reaction rate constants
global k_minus
global k_plus
global k1
global k2
global n

%assign parameter values
% Production and degradation terms
alpha_mrnaC=1.0;
alpha_mrnaY=1.0;
alpha_Rg=1.0;
alpha0_mrnaY=1.0;
beta_C=1.0;
gamma_C=1.0;
gamma_mrnaC=1.0;
gamma_Rg=1.0;
gamma_mrnaY=1.0;

% Reaction rate constants
k_minus=1.0;
k_plus=1.0;
k1=1.0;
k2=1.0;
n=2.5;

% Set simulation parameters
ODEFUN=@crispr_ode;
options=odeset('Refine', 6);
Tend=300;

% Set initial condition: state = [mrnaC, C, Rg, mrnaY]
x0=[0.2, 0.1, 0.3, 0.1]';

% Run simulation
[t,S]=ode45(ODEFUN, [0,Tend], x0, options);

% Make a plot
figure(1)
hold on
plot(t, S(:,1), 'c', 'Linewidth', 1.5)
plot(t, S(:,2), 'b', 'Linewidth', 1.5)
plot(t, S(:,3), 'r', 'Linewidth', 1.5)
plot(t, S(:,4), 'k--', 'Linewidth', 1.5)
legend('mrnaC', 'C', 'Rg', 'mrnaY', 'Location','Best')
xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')

end

% ODE model
function dS=crispr_ode(t,s)

    % Production and degradation terms
    global alpha_mrnaC
    global alpha_mrnaY
    global alpha_Rg
    global alpha0_mrnaY
    global beta_C
    global gamma_C
    global gamma_mrnaC
    global gamma_Rg
    global gamma_mrnaY

    % Reaction rate constants
    global k_minus
    global k_plus
    global k1
    global k2
    global n

    % Differential terms
    mrnaC=s(1); % dCas9 mRNA
    C=s(2); % dCas9
    Rg=s(3); % sgRNA
    mrnaY=s(4); % mRNA to produce YFP

    % Non-differential terms from QSSA
    K = k_plus/k_minus;
    K_a = k2/k1;
    CRg = K*C*Rg;

    mrnaC_dt= alpha_mrnaC - gamma_mrnaC*mrnaC;
    C_dt= beta_C*mrnaC - gamma_C*C + k_minus*CRg - k_plus*C*Rg;
    Rg_dt=alpha_Rg - gamma_Rg*Rg + k_minus*CRg - k_plus*C*Rg;
    mrnaY_dt=alpha_mrnaY/(1 + (CRg/K_a)^n) - gamma_mrnaY*mrnaY + alpha0_mrnaY;

    dS=[mrnaC_dt, C_dt, Rg_dt, mrnaY_dt]';
end

