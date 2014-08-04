% Not sure if this makes sense- looking for feedback!

% Very tiny script to fit the three SarA P1 data points from Figure 2
% in Malone et al., http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2693297/
% Disregarding the timepoint at 16h because it's obviously reached
% steady-state at that point (which also, clearly, means that a linear
% approximation is incorrect, but it allows us to define a single value for
% SarA P1 expression)

% Values for fluorescence were measured super-hackily by interpolating the
% pixels in the figure (which is also why they are to the nearest 50)

time_h=[4 5 6];
fluorescence_rel=[1850 6000 9850]; % (measured relative to cell density)
fluorescence_rel_rel = fluorescence_rel ./ fluorescence_rel(1);

stats=regstats(fluorescence_rel_rel,time_h,'linear')
coefficients=stats.beta