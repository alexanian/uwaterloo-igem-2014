% Initialize
clc

dx = 0.01; % Width of each point in cm
dy = 0.01; % Height of each point in cm
tend = 24; % Time to stop the simulation in hours

% Evaluate
matlab_filename = 'PDE_implementation';

expression = ['[D,R,T,N,dt,times] = ',matlab_filename,'(dx,dy,tend);'];
eval(expression);

% Plot
num_steps = length(N(1,1,:));
dim = 5;
sample_rate = round(num_steps/(dim-1))-1

C = max(max(max(N))); % used to get the scales consistent

% Main loop for plot
for num = 1:dim
    
    if num == 1
        nn = 1;
    else
        nn = (num-1)*sample_rate+1;
    end
    
    a = nn*dt;
    
    figure(1)
    
    % Top row shows Donors
    subplot(3,dim,num)
    pcolor(D(:,:,nn))
    caxis([0 C])
    title(strcat('Donors, t=',num2str(a),' hours'))
    shading interp
    
    % Middle row shows Recipients
    subplot(3,dim,num+dim)
    pcolor(R(:,:,nn))
    caxis([0 C])
    title(strcat('Recipients, t=',num2str(a),' hours'))
    shading interp
    
    % Third row shows Transconjugants
    subplot(3,dim,num+2*dim)
    pcolor(T(:,:,nn))
    caxis([0 C])
    title(strcat('Transconjugants, t=',num2str(a),' hours'))
    shading interp
    
end

total_rec = zeros(1,num_steps);
total_cel = zeros(1,num_steps);
total_other = zeros(1,num_steps);
frac_rec  = zeros(1,num_steps);

for nn = 1:num_steps
    total_rec(nn) = sum(sum(R(:,:,nn)));
    total_cel(nn) = sum(sum(N(:,:,nn)));
    total_other(nn) = sum(sum(D(:,:,nn)+T(:,:,nn)));
    frac_rec(nn) = total_rec(nn)/total_cel(nn);
end

scale = (1:num_steps)*dt;
length(scale);
length(total_other);

figure(2)
% subplot(2,1,1)
plot(scale,total_other,scale,total_rec,'linewidth',1.5)
title('Cell populations over time (1cm by 1cm grid)')
legend('Donors/Transconjugants','Recipients')
axis([0 num_steps*dt 0 1.15*max(total_cel)])
xlabel('Time (h)')
ylabel('Number of cells')
