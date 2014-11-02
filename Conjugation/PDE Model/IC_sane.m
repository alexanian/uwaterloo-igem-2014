% This function determines the Initial conditions of the system

function P = IC_sane(xlen,ylen,K)
% xlen and ylen specify number of points in x or y direction
% K is carrying capacity

init_grid = zeros(xlen,ylen);   % Initialization

% Loop through each point on the grid
for ii = 1:xlen
    for jj = 1:ylen
        r = rand; % Add in some randomness to make things interesting
        init_grid(ii,jj) = K*r; % The distribution
    end
end

P = init_grid;

end