function [t,R] = sRNA_model_mechanism_1
	% Species
    titles = {'sRNA'; 'mRNA'; 'YFP'; 'Hfq-mRNA'; 'Hfq'; 'Hfq-sRNA Complex'; 'Hfq-mRNA-sRNA Complex'};
	nEqns = 7;
    
	% Initial Conditions
    % ** Note: These were obtained from running a simulation to steady
    %          state with no production of sRNA within it. That's why they
    %          seem oddly specific.
	R0=zeros(nEqns,1);
    R0(1)=0; % sRNA
    R0(2)=0.721347520829037; % mRNA
    R0(3)=18.6892224137109; % YFP
    R0(4)=0.721347520829037; % Hfq-mRNA
    R0(5)=10070.4043649796; % Hfq
    R0(6)=0; % Hfq-sRNA Complex
	R0(7)=0; % Hfq-sRNA-mRNA Complex
	Tend = 2*24*60*60; % Time
	
	% ODE Solver and Plots
	%opt = odeset('InitialStep',0.001);
	%[t, R] = ode45(@packddt, [0,Tend], R0); % opt);
    [t, R] = ode15s(@packddt, [0,Tend], R0); % opt);
    
	figure(1);
    for i = 1:length(titles)
        subplot(2, 4, i);
        plot(t/3600, R(:,i), 'LineWidth', 2);
        title(titles(i));
        axis tight
    end
    hold off;
    

    
end

function dR = packddt(~, R)
	
    % Constants 
    % * rate of productions - alpha
    % * rate of degradation - beta
	% * lowercase means RNA, uppercase mean protien product
    alpha_h = 1/600; % nM/s
    beta_h = 0.002310490602; % nM/s
    alpha_H = 0.9; % translation not provided?
    beta_H = 6.42*10^(-5);  % degradation not provided?
    
    alpha_s = 1/600; % nM/s
    beta_s = 0.002310490602; % nM/s
    
    alpha_m = 1/600;
    beta_m = 0.002310490602;
    alpha_M = 1/600;
    beta_M = 6.42*10^(-5);
    
    
    k1 = 10^(6);
    k_1 = 0.7*10^(-4); % the backward rate from Hfq-sRNA association was taken 
    %             (previously) to be zero. Hence this parameter didn't
    %             exist.
    k2 = 3.5*10^(6);
    
    % degradosome is taken to be mass action
    % ** NOTE: This is where the major difference is.
    %          In a gram-negative cell there is a specialized 
    %          unit called a degradosome whose job it is to break
    %          down the Hfq-sRNA-mRNA complex. In gram-negative there isn't
    %          one, so the degradation is taken to be simply mass action.
    %          Without a degradosome, and with no data regarding the
    %          degradation of the Hfq-sRNA-mRNA complex, I'm not entirely
    %          sure what to do here.
    %    TL;DR: This should be mass action, and not Michaelis-Menton. This
    %           is where the confusion was.
    %V_DSM = 1;
    %Km_DSM = 1;
    k3 = 0.7*10^(-4);
    
	% Assign biological names to elements of R vector
    s = R(1);   % sRNA
    m = R(2);   % YFP-mRNA
    M = R(3);   % YFP
    h = R(4);   % Hfq mRNA
    H = R(5);   % Hfq
    Hs = R(6);  % Hfq-sRNA Complex
    Hms = R(7); % Hfq-mRNA-sRNA Complex
    
    
    % sRNA ODE's
    d_s = alpha_s - beta_s*s - k1*H*s + k_1*Hs;
    d_m = alpha_m - beta_m*m - k2*Hs*m;
    d_M = alpha_M*m - beta_M*M;
    d_h = alpha_h - beta_h*h;
    d_H = alpha_H*h - beta_H*H + k3*Hms;% (V_DSM*Hms)/(Hms + Km_DSM) - k1*H*s;
    d_Hs = k1*H*s - k2*Hs*m - k_1*Hs;
    d_Hms = k2*Hs*m - k3*Hms; % (V_DSM*Hms)/(Hms + Km_DSM);
    
	
	dR=[d_s; d_m; d_M; d_h; d_H; d_Hs; d_Hms];
	
end



















