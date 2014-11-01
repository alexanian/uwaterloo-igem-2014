function Sout = RNAi_System_YFP_Output( Sin )
% Extracts repressed YFP statistic from CRISPRi_Model output
    Sout = Sin(:,3);
end