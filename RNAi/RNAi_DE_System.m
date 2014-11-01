function dS = RNAi_DE_System( ~, S, p )
% RNAi_DE_System del of RNAi system to be used with ode45/ode23
    % Parameters & Explanation
    % * rate of productions - alpha
    % * rate of degradation - beta
	% * lowercase means RNA, uppercase mean protien product
    alpha_h = p(1); % nM/s
    beta_h = p(2); % nM/s
    alpha_H = p(3); % translation not provided?
    beta_H = p(4);  % degradation not provided?
    alpha_s = p(5); % nM/s
    beta_s = p(6); % nM/s 
    alpha_m = p(7);
    beta_m = p(8);
    alpha_M = p(9);
    beta_M = p(10);
    k1 = p(11);
    k_1 = p(12);
    k2 = p(13);
    k3 = p(14);

    % Differential terms from elements of R vector
    s = S(1);   % sRNA
    m = S(2);   % YFP-mRNA
    M = S(3);   % YFP
    h = S(4);   % Hfq mRNA
    H = S(5);   % Hfq
    Hs = S(6);  % Hfq-sRNA Complex
    Hms = S(7); % Hfq-mRNA-sRNA Complex
    
    if S(6)>0.06
        8;
    end
    
    % RNAi Differential System
    d_s = alpha_s - beta_s*s - k1*H*s + k_1*Hs;
    d_m = alpha_m - beta_m*m - k2*Hs*m;
    d_M = alpha_M*m - beta_M*M;
    d_h = alpha_h - beta_h*h;
    d_H = alpha_H*h - beta_H*H + k3*Hms+ k_1*Hs; % (V_DSM*Hms)/(Hms + Km_DSM) - k1*H*s;
    d_Hs = k1*H*s - k2*Hs*m - k_1*Hs - beta_H*Hs;
    d_Hms = k2*Hs*m - k3*Hms - beta_H*Hms;    % (V_DSM*Hms)/(Hms + Km_DSM);
    
	%%added terms Hs -> H, s
    %%degradation term of Hs (assuming when it meets a protieosome it chews
    %%up the whole complex using the rate beta_H that was representative of
    %%the protieosome
    %%degradation term of Hms using the same assumptions as above
    
	dS=[d_s; d_m; d_M; d_h; d_H; d_Hs; d_Hms];
end

