function Parabolic
    % Chao Qian
    % first initialize some parameters  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the advection speed
    a=1.;
    %
    ete = 2.;
    % number of grid points in space (not counting the point at x=1, which is the same as the point at x=0, due to the periodicity)
    N=80;
    % spatial step
    h =(20./(N));
    % two relevant time step size
    dt_adv = h/a;
    dt_diff = h^2/(2*ete);
    % timestep limit 
    dt1 = 0.9 * 0.5 * min(dt_adv,dt_diff);
    dt2 = 0.9*dt_adv;
    % time at which we want to end the simulation
    t_end=1.;
    % number of timesteps to be taken
    it1=t_end/dt1;
    it2=t_end/dt2;
    % number of different methods we want to try
    n_methods=2;
    
    % initialize some arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spatial grid (no need to include the boundary, since the solution remains 0.)
    x=(-10:h:10-h);
    v_new=zeros(N,n_methods);
    v_old=zeros(N,n_methods);
    % the initial condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:n_methods
       v_old(1,i) = 0;
       
       v_old(2:N,i)=1/sqrt(0.01)*exp(-1 * x(2:N).^2/(4*ete*0.01));
    end
    
    % get the matrices for the implicit method
    % solve for A V(n+1) = B V(n)
    [A,B]=semi_implicit(N,a,h,dt2,ete);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iter = 1:it1
        % method 1: explicit method
        
        v_new(2:N-1,1)=v_old(2:N-1,1)-a*dt1/h*(v_old(2:N-1,1)-v_old(1:N-2,1))+ete*dt1/h^2*(v_old(3:N,1)-2*v_old(2:N-1,1)+v_old(1:N-2,1));
        v_new(1,1)=0;
        v_new(N,1) = v_old(N,1)-a*dt1/h*(v_old(N,1)-v_old(N-1,1))+ete*dt1/h^2*(v_old(1,1)-2*v_old(N,1)+v_old(N-1,1));
         
        v_old(:,1) = v_new(:,1);
    end
    
    for iter = 1:it2
        % method 2: soemi-implicit method
        v_new(:,2)=A\(B*v_old(:,2));
        v_old(:,2) = v_new(:,2);
    end
    % graphical output 
    plot(x,v_new(:,1),'*k-')
    hold on 
    
    plot(x,v_new(:,2),'>b-')
    xlabel('x')
    ylabel('v')
    title('Parabolic')
    legend('Emplicit','Semi-implicit','Location', 'NorthWest','FontSize',8);
    hold off
   
end 


function [C,D]=semi_implicit(N,a,h,dt2,ete)
  e = ones(N,1);
  C = spdiags([-ete * dt2/h^2*e (1+ 2 * ete * dt2 /h^2)*e -ete * dt2/h^2*e], [-1 0 1], N, N);
  C(N,1)=-ete*dt2/h^2;
  C(1,2) = 0;
  C(1,1) = 1;
  
  D= spdiags([a*dt2/h*e (1-a*dt2/h)*e],[-1 0],N,N);
  D(1,1) = 1;
end