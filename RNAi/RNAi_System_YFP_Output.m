function Sout = RNAi_System_YFP_Output( S_in )
% Extracts repressed YFP statistic from RNAi_DE_System output
    Sout = S_in(:,3);
end