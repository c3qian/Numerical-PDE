% Initialization of table variables
execut_time1 = zeros(9,1);
execut_time2 = zeros(9,1);
timestep1 = zeros(9,1);
timestep2 = zeros(9,1);
timesize1 = zeros(9,1);
timesize2 = zeros(9,1);
advect_tsz = zeros(9,1);
diffus_tsz = zeros(9,1);

N = [11,21,41,81,161,321,641, 1281, 2561]';

for i = 1:9
    [timesize1(i),timesize2(i),advect_tsz(i),diffus_tsz(i),execut_time1(i),execut_time2(i),timestep1(i),timestep2(i)] = advdiff(2^(i-1)*10);
end

my_table = table(N, execut_time1,  timestep1, timesize1, execut_time2, timestep2, timesize2,  advect_tsz, diffus_tsz);

function  [dt1,dt2,dt_adv,dt_diff,exe_t1,exe_t2,n_it1,n_it2]= advdiff(N)
    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the advection speed
    a=1.;
    %
    eta = 2.;
    % spatial step
    h =(20./(N));
    % two relevant time step size
    dt_adv = h/a;
    dt_diff = h^2/(2*eta);
    % timestep limit 
    dt1 = 0.9*0.5*min(dt_adv,dt_diff);
    dt2 = 0.9*dt_adv;
    % time at which we want to end the simulation
    t_end=1.;
    % number of timesteps to be taken
    n_it1=(t_end-0.01)/dt1;
    n_it2=(t_end-0.01)/dt2;
    % number of different methods we want to try
    n_methods=2;
    
    % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spatial grid (no need to include the boundary, since the solution remains 0.)
    x=(-10:h:10-h);
    % arrays for the numerical approximations at times n and n+1
    v_new=zeros(N,n_methods);
    v_old=zeros(N,n_methods);
    % the initial condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:n_methods
       v_old(1,i) = 0;
       v_old(2:N,i)=1/sqrt(0.01)*exp(-x(2:N).^2/(4*eta*0.01));
    end
    
    % get the matrices for the implicit method
    % A V(n+1) = B V(n)
    [A,B]=mat_advdiff_si(N,a,h,dt2,eta);
    
    tic
    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    for iter = 1:n_it1
        % method 1: explicit method
        v_new(2:N-1,1)=v_old(2:N-1,1)-a*dt1/h*(v_old(2:N-1,1)-v_old(1:N-2,1))+eta*dt1/h^2*(v_old(3:N,1)-2*v_old(2:N-1,1)+v_old(1:N-2,1));
        v_new(1,1)=0;
        v_new(N,1) = v_old(N,1)-a*dt1/h*(v_old(N,1)-v_old(N-1,1))+eta*dt1/h^2*(v_old(1,1)-2*v_old(N,1)+v_old(N-1,1));
         
        v_old(:,1) = v_new(:,1);
    end
    exe_t1 = toc;
    tic
    for iter = 1:n_it2
        % method 2: semi-implicit method
        v_new(:,2)=A\(B*v_old(:,2));
        v_old(:,2) = v_new(:,2);
    end
    exe_t2 = toc;
end 
% the matrices for the implicit method
%-------------------------
% implicit central in space (BC)
function [A,B]=mat_advdiff_si(N,a,h,dt2,eta)
  e = ones(N,1);
  A = spdiags([-eta*dt2/h^2*e (1+2*eta*dt2/h^2)*e -eta*dt2/h^2*e], [-1 0 1], N, N);
  A(1,1) = 1;
  A(1,2) = 0;
  A(N,1)=-eta*dt2/h^2;
  B= spdiags([a*dt2/h*e (1-a*dt2/h)*e],[-1 0],N,N);
  B(1,1) = 1;
end