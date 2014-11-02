function dS = RNAi_DE_System( ~, S, p )
% RNAi_DE_System of RNAi system to be used with ode23s/ode15s
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
    k_1 = p(11);
    k_m1 = p(12);
    k_2 = p(13);
    k_m3 = p(14);
    k_4 = p(15);
    beta_ms = p(16);

    % Differential terms from elements of R vector
    s = S(1);   % sRNA
    m = S(2);   % YFP-mRNA
    M = S(3);   % YFP
    h = S(4);   % Hfq mRNA
    H = S(5);   % Hfq
    Hs = S(6);  % Hfq-sRNA Complex
    Hms = S(7); % Hfq-mRNA-sRNA Complex
    ms = S(8);  % bound (inactive) mRNA-sRNA

    if S(6)>0.06
        8;
    end
    
    % RNAi Differential System
    %      production      degradation     Hs complex         Hms complex   ms complex
    d_s =  alpha_s         -beta_s*s       -k_1*H*s +k_m1*Hs               -k_4*m*s;
    d_m =  alpha_m         -beta_m*m                          -k_2*Hs*m    -k_4*m*s;
    d_M =  alpha_M*m       -beta_M*M;
    d_h =  alpha_h         -beta_h*h;
    d_H =  alpha_H*h       -beta_H*H       -k_1*H*s +k_m1*Hs  +k_m3*Hms;
    d_Hs =                 -beta_H*Hs      +k_1*H*s -k_m1*Hs  -k_2*Hs*m;
    d_Hms=                 -beta_H*Hms                        +k_2*Hs*m -k_m3*Hms;
    d_ms = k_4*m*s         -beta_ms*ms                        -k_m3*Hms;
    
	dS=[d_s; d_m; d_M; d_h; d_H; d_Hs; d_Hms; d_ms]; 
end

