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
alpha_mrnaC=0.05;%0.0011;
alpha_mrnaY=0.05;%0.0011;
alpha_Rg=0.01;%0.0011;
beta_C=(0.0057 + 0.4797)/2; % midpoint of range
gamma_C=-5.6408*10^-4;
gamma_B=-5.6408*10^-4;
gamma_mrnaC=0.1734;
gamma_Rg=0.1734;
gamma_mrnaY=0.1734;

% Reaction rate constants
k_minus=0.1; % midpoint of range
k_plus=0.1; % midpoint of range
K_a=0.28;
n=2.5;

% Set simulation parameters
ODEFUN=@CRISPRi_System;
options=odeset('Refine', 6);
Tend=300;
params = [alpha_mrnaC alpha_mrnaY alpha_Rg beta_C...
    gamma_C gamma_B gamma_mrnaC gamma_Rg gamma_mrnaY...
    k_minus k_plus K_a n];

% Set initial condition: state = [mrnaC, C, Rg, mrnaY]
x0=[0, 0, 0, 0.29, 0.29]';

% Run simulation
[t,S]=ode45(ODEFUN, [0,Tend], x0, options, params);

% Algebraic term from QSSA
B = k_plus*S(:,2).*S(:,3)/(k_minus+gamma_B);

% Make a plot
figure;
hold on
plot(t, S(:,1), 'color', [0.77 0.85 0.89], 'Linewidth', 1.5)
plot(t, S(:,2), 'color', [0.39 0.58 0.93], 'Linewidth', 1.5)
plot(t, S(:,3), 'color', [0.88 0.24 0.19], 'Linewidth', 1.5)
plot(t, B, 'color', [0.42 0.19 0.51], 'Linewidth', 1.5)
plot(t, S(:,4), 'color', [0.98 0.93 0.36], 'Linewidth', 1.5)
legend('dCas9 mRNA', 'dCas9', 'sgRNA', 'dCas9-sgRNA Complex',...
    'YFP mRNA', 'Location','Northeast')
xlabel('Time (minutes)')
ylabel('Concentration (nM)')
ylim([0 0.5])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.001 .001] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:0.5, ...
  'LineWidth'   , 1         );

figure;
hold on
plot(t, S(:,1), 'color', [0.77 0.85 0.89], 'Linewidth', 1.5)
plot(t, S(:,2), 'color', [0.39 0.58 0.93], 'Linewidth', 1.5)
plot(t, S(:,3), 'color', [0.88 0.24 0.19], 'Linewidth', 1.5)
plot(t, B, 'color', [0.42 0.19 0.51], 'Linewidth', 1.5)
plot(t, S(:,5), 'color', [0.98 0.93 0.36], 'Linewidth', 1.5)
legend('dCas9 mRNA', 'dCas9', 'sgRNA', 'dCas9-sgRNA Complex',...
    'YFP mRNA', 'Location','Northeast')
xlabel('Time (minutes)')
ylabel('Concentration (nM)')
ylim([0 0.5])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.001 .001] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:0.5, ...
  'LineWidth'   , 1         );
end