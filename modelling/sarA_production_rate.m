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

% Then normalized the observed change to get an accurate rate
fluor_norm = (fluor_molec - fluor_molec(1)) / fluor_molec(1);

linearCoef = polyfit(time_min,fluor_norm,1);
linearFit = polyval(linearCoef,time_min);
plot(time_min, fluor_norm, 'k.', time_min,linearFit, 'r');
xlabel('Time (min)'); ylabel('Normalized change in # of Fluorescent Molecules');
xlim([180 600]); ylim([0 0.6]);

sarA_per_molecule_rate = linearCoef(1)
