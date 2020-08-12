function Goda 
    v = 16;
    a=1.;
    % number of grid points in space (not counting the point at x=1, which is the same as the point at x=0, due to the periodicity)
    N=50;
    % spatial step
    h =(1./(N));
    % safety constant for stability (should be smaller than 1)
    cfl=.8;
    % CFl timestep limit for the explicit methods
    dt=cfl* h/a;
    % time at which we want to end the simulation
    t_end=1.;
    % number of timesteps to be taken
    n_it = t_end/dt;
    % number of different methods we want to try
    n_methods=1;
    
    
        % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spatial grid (not including the point at x=1, which is the same as the point at x=0, due to the periodicity)
    x=(0:h:1-h);
    % temporal grid
    t=(0:dt:n_it*dt);
    % arrays for the numerical approximations at times n and n+1
    v_new=zeros(N,n_methods);
    v_old=zeros(N,n_methods);
    %v_new_lw=zeros(N,n_methods);
    
    % array for the exact solution
    v_exact=zeros(N,1);

    % the initial condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:n_methods
       %problem 1:
       %v_old(:,i)=exp(-50*(x-.5).^2);
       v_old(:,i)= -8/(cosh(2*x).^2);
    end
    
    [C,D]=mat_linadv_CN(N,a,h,dt);
    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:n_it
        v_new(:,2)=A\(B*v_old(:,2));
        v_exact(:) = -0.5*v./(cosh(0.5*sqrt(v)*(x-v*t(iter+1)-x_0))).^2;
              
        plot(x,v_exact(:),'^r-')
        hold on
        
        for i=1:n_methods
           plot(x,v_new(:,2),'+g-')
           
        end
        axis([0 1 -.2 1.2])
        xlabel('x')
        ylabel('v')
        title('linear advection simulation')
        hold off
        pause(0.001)
		
        % prepare for the next iteration
        v_old=v_new;
    end        
end


function [C,D]=mat_linadv_CN(N,a,h,dt)
  e = ones(N,1);
  para1 = -0.5*dt/h^3*e;
  para2 = dt/h^3*e-dt/h(v);
  para3 = -para2;
  para4 = 0.5*dt/h^3*e;
  C = spdiags([para1 para2 e para3 para4], [-2 -1 0 1 2], N, N);
  %C(1,N)=-0.25*a*dt/h;
  %C(N,1)=0.25*a*dt/h;
  D = speye(N);
end