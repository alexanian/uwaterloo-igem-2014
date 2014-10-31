function [D,R,T,N,dt,times] = PDE_implementation(dx,dy,tend)

format long

xlen = 1; % Width of the simulated area in cm
ylen = 1; % Height of the simulated area in cm
inmodx = dx/xlen; % Inverse number of xpoints
inmody = dy/ylen; % Inverse number of ypoints

d  = 0.001;
r  = 0.1;

dt = r*(inmodx^2+inmody^2)/d; % Size of time step in hours

doub = 0.5;%0.25; % Doubling time in hours
growth = (0.69314718056)/doub; % Growth rate
rN = growth;

cells_per_square_meter = 2.8867513459 * 10^11;
K = cells_per_square_meter*dx*dy/10000 % total number of cells in a point (10micronx10micron square)
lambda = 0.2/K;%(10^(-4))/K

dim = 5;
x_total = round(xlen/dx);
y_total = round(ylen/dy);
num_steps = round(tend/dt)
sample = round(num_steps/(dim-1))-1;

D = zeros(x_total,y_total,num_steps);
R = zeros(x_total,y_total,num_steps);
T = zeros(x_total,y_total,num_steps);
N = zeros(x_total,y_total,num_steps);

R(:,:,1) = IC_sane(x_total,y_total,K)/2;
D(:,:,1) = IC_sane(x_total,y_total,K)/2;

times = zeros(1,num_steps);

num = 0;
tic
for nn = 1:num_steps-1
    N(:,:,nn) = D(:,:,nn)+R(:,:,nn)+T(:,:,nn);
    for ii = 1:x_total
        for jj = 1:y_total
            if ii == 1
                iL = x_total;
                iR = 2;
            elseif ii == x_total;
                iL = x_total-1;
                iR = 1;
            else
                iL = ii-1;
                iR = ii+1;
            end
            
            if jj == 1
                jL = x_total;
                jR = 2;
            elseif jj == x_total;
                jL = x_total-1;
                jR = 1;
            else
                jL = jj-1;
                jR = jj+1;
            end
            
            dN2 = ((N(iR,jj,nn)-2*N(ii,jj,nn)+N(iL,jj,nn))/(dx^2)+(N(ii,jR,nn)-2*N(ii,jj,nn)+N(ii,jL,nn))/(dy^2))*d*dt/(N(ii,jj,nn));
            
            gro = (1-N(ii,jj,nn)/K)*dt;
            conj = (D(ii,jj,nn)+T(ii,jj,nn))*R(ii,jj,nn)*dt*lambda;
            multer = 1+dN2+rN*gro;
            
            D(ii,jj,nn+1) = D(ii,jj,nn)*multer;
            R(ii,jj,nn+1) = R(ii,jj,nn)*multer-conj;
            T(ii,jj,nn+1) = T(ii,jj,nn)*multer+conj;
        end
    end
    
    if nn == num*sample+1
        num = num+1
    end
    
    times(nn+1) = times(nn)+dt;
    
end

times

toc

end