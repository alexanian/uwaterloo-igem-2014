% Estimation of mRNA production from SarA promoter
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

% Curve fitting of exponent
ft=fittype('exp1');
expCoeff = coeffvalues(fit(time_min',fluor_molec',ft));
expFit = expCoeff(1)*exp(expCoeff(2)*time_min);

h = plot(time_min, fluor_molec, 'kd', time_min,expFit, 'b');
set(h, {'MarkerFaceColor'}, {'k'; 'b';}, {'MarkerSize'}, {6; 1;},...
    {'LineWidth'},{1; 2;})
xlabel('Time (min)'); ylabel('# of Fluorescent Molecules');
xlim([180 600]); 

sarA_per_molecule_rate = expCoeff(2)

% Estimation of dCas9 protein degradation rate
% Half-life of SarA protein in S Aureus given in Michalik et al
% http://www.mcponline.org/content/11/9/558.abstract
half_life_h = 20.48;
gamma_C = -log(2)/(half_life_h*60);

% Estimation of dCas9 translation rate
peptide_s = [0.59 3.17];
dCas9_s = peptide_s / 1372;
beta_C = dCas9_s * 60;