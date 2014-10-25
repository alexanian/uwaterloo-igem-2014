% Code describing a number of estimates made from the literature- all
% resulting parameters are found in CRISPRi_Literature_Parameters.mat
%
% Estimation of mRNA production from SarA promoter
%%%
% Fluorescence (OD) read from Cheung et al, see figure 4B in
% http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2258893/
% Values for relative fluorescence measured super-hackily by interpolating the
% pixels in the figure (which is also why they are to the nearest 10)
time_h = [3 4 5 6 7 8 9 10];
time_min = time_h * 60;
fluor_OD_sarAneg = [11680 13260 12480 12400 13670 15860 17120 18290];
fluor_OD_sarA = [9500 10590 10700 11310 12420 13400 14000 15080];

% I then converted from fluorescence to number of molecules
% using Wu & Pollard: 0.00676 A.U / mean # molecules per cell
fluor_molec = fluor_OD_sarA / 0.00676;
fluor_nanomol = fluor_molec / (6.0331415e23 * 10e-9);

% Convert to a per-concentration rate
staph_cell_volume_l = 4/3*pi*1.1e-6^3 * 1000;
fluor_nM = fluor_nanomol / staph_cell_volume_l;

% Curve fitting of exponent
ft=fittype('exp1');
expCoeff = coeffvalues(fit(time_min',fluor_nanomol',ft));
expFit = expCoeff(1)*exp(expCoeff(2)*time_min);
sarA_per_mol_rate = expCoeff(2);

figure;
h = plot(time_min, fluor_nanomol, 'kd', time_min,expFit, 'b');
set(h, {'MarkerFaceColor'}, {'k'; 'b';}, {'MarkerSize'}, {8; 1;},...
    {'LineWidth'},{1; 2;})
xlabel('Time (min)'); ylabel('Fluorescent Molecules in Culture (nanomoles)');
xlim([180 600]); 

%%%
% Estimation of dCas9 protein degradation rate
%%%
% Half-life of SarA protein in S Aureus given in Michalik et al
% http://www.mcponline.org/content/11/9/558.abstract
halfLife_h = 20.48;
gamma_C = -log(2)/(halfLife_h*60);

%%%
% Estimation of dCas9 translation rate
%%%
peptide_s = [0.59 3.17];
dCas9_s = peptide_s / 1372;
dCas9_m = dCas9_s * 60;
beta_C = [dCas9_m(1)*0.22 dCas9_m(2)*3.46]