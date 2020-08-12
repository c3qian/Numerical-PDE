    v = 16;
    x0 = 0;    
    dx = 0.1;
    h = dx;
    N = 16 / dx;
    u_max = -8;
    dt = dx / ((4 / dx ^ 2) + 6 * abs(u_max)) * 0.5;

    t_end = 1;
    n_it = t_end / dt;
    
    % number of different methods we want to try
    n_methods = 2;   
    v_new = zeros(N, n_methods);
    temp = zeros(N, n_methods); 
    v_old = zeros(N, n_methods);
    dia1 = zeros(N, n_methods);

    x = (- 8 : dx : 8 - dx);
    t=(0 : dt : n_it*dt);
    v_exact = zeros(N, 1);
    
    % leapfrog method: 
    % initialization
    for i = 1 : n_methods
       %P b
       v_old(:,i)=-8./cosh(2*x).^2;
       %P c
       %v_old(:,i)=-8./cosh(2*x+8).^2-2./cosh(x).^2;
       %P d
       %v_old(:,i)=-5* exp(-(x/2).^2);
    end
    
    
    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the first iteration    
  
    v_new(3:N-2, 1) = v_old(3:N-2, 1) + 2 * dt * ((v_old(2:N-3, 1) + v_old(3:N-2, 1) + ...
            +v_old(4:N-1, 1)) .* (v_old(4:N-1, 1) - v_old(2:N-3, 1)) / (dx) - (v_old(5:N, 1) - 2 * v_old(4:N-1, 1) +...
            + 2 * v_old(2:N-3, 1) - v_old(1:N-4, 1)) / (2 * dx ^ 3) );

    
    v_new(1, 1) = v_old(1, 1) + dt * ((v_old(N, 1) + v_old(1, 1) + ...
            + v_old(2, 1)) * (v_old(2, 1) - v_old(N, 1)) / (dx) -(v_old(3, 1) - 2 * v_old(2, 1) +...
            + 2 * v_old(N, 1) - v_old(N - 1, 1)) / (2 * dx ^ 3) );
    v_new(2, 1) = v_old(2, 1) + dt * ((v_old(1, 1) + v_old(2, 1) + ...
            + v_old(3, 1)) * (v_old(3, 1) - v_old(1, 1)) / (dx) -(v_old(4, 1) +...
            - 2 * v_old(3, 1) + 2 * v_old(1, 1) - v_old(N, 1)) / (2 * dx ^ 3) );
    v_new(N - 1, 1) = v_old(N - 1, 1) + dt * ((v_old(N - 2, 1) + ...
            +v_old(N - 1, 1) + v_old(N, 1)) * (v_old(N, 1) - v_old(N - 2, 1)) / (dx) + ...
            -(v_old(1, 1) - 2 * v_old(N, 1) + 2 * v_old(N - 2, 1) - v_old(N - 3, 1)) / (2 * dx ^ 3) );
    v_new(N, 1) = v_old(N, 1) + dt * ((v_old(N - 1, 1) + v_old(N, 1) +...
            + v_old(1, 1)) * (v_old(1, 1) - v_old(N - 1, 1)) / (dx) -(v_old(2, 1) +...
            - 2 * v_old(1, 1) + 2 * v_old(N - 1, 1) - v_old(N - 2, 1)) / (2 * dx ^ 3) );
    v_olld =v_old(:,1);
       
    % after the first iteration
    for iter = 1: n_it 
        % method 1: leapfrog scheme
        
        v_new(3:N-2, 1) = v_olld(3:N-2, 1) + 2 * dt * ((temp(2:N-3, 1) + temp(3:N-2, 1) + ...
            +temp(4:N-1, 1)) .* (temp(4:N-1, 1) - temp(2:N-3, 1)) / (dx) - (temp(5:N, 1) - 2 * temp(4:N-1, 1) +...
            + 2 * temp(2:N-3, 1) - temp(1:N-4, 1)) / (2 * dx ^ 3) );      
        v_new(1, 1) = v_olld(1, 1) + 2 * dt * ((temp(N, 1) + temp(1, 1) + temp(2, 1)) * (temp(2, 1) - temp(N, 1)) / (dx - (temp(3, 1) - 2 * temp(2, 1) + 2 * temp(N, 1) - temp(N - 1, 1)) / (2 * dx ^ 3) ));
        v_new(2, 1) = v_olld(2, 1) + 2 * dt * ((temp(1, 1) + temp(2, 1) + temp(3, 1)) * (temp(3, 1) - temp(1, 1)) / (dx)+...
            - (temp(4, 1) - 2 * temp(3, 1) + 2 * temp(1, 1) - temp(N, 1)) / (2 * dx ^ 3) );
        v_new(N - 1, 1) = v_olld(N - 1, 1) + 2 * dt * ((temp(N - 2, 1) + temp(N - 1, 1) + temp(N, 1)) * (temp(N, 1) + ...
            - temp(N - 2, 1)) / (dx) - (temp(1, 1) - 2 * temp(N, 1) + 2 * temp(N - 2, 1) - temp(N - 3, 1)) / (2 * dx ^ 3) );
        v_new(N, 1) = v_olld(N, 1) + 2 * dt * ((temp(N - 1, 1) + temp(N, 1) + temp(1, 1)) * (temp(1, 1) - temp(N - 1, 1)) / (dx)+ ...
            - (temp(2, 1) - 2 * temp(1, 1) + 2 * temp(N - 1, 1) - temp(N - 2, 1)) / (2 * dx ^ 3) );
        
        % method 2: Goda scheme
        
        [C,D]= Goda(N,v_old(:,2),h,dt);    
        v_new(:,2)=C\(D*v_old(:,2));
        v_old = v_new;
        
        % the exact solution with        
        v_exact(:)=-8./cosh(2*(x-16*iter*dt)).^2;
        
        
        % graphical output      
        plot(x, v_exact(:),'^r-')         
        hold on
        plot(x,v_exact,'*k-')
        plot(x,v_old(:,1),'>b-')
        plot(x,v_old(:,2),'+r-')
        xlabel('x')
        ylabel('v')
        title('t=1 Gaussian profile')
        legend('exact','LF','Goda','location','SouthEast','FontSize',5)
  
        hold off            

    end
    
function [C,D]=Goda(N,v_old, h,dt)
   e = ones(N,1);
  a=[v_old(N)+v_old(1);v_old(1:N-1)+v_old(2:N)];
  b=[v_old(2:N)+v_old(1:N-1);v_old(1)+v_old(N)];
  C = spdiags([-dt/2/h^3*e dt/h^3+dt/h*b e -dt/h^3-dt/h*a dt/2/h^3*e], [-2 -1 0 1 2], N, N);
  C(1,N-1)=- dt/2/h^3;
  C(1,N)= dt/h^3+dt/h*b(N);
  C(2,N)= -dt/2/h^3;
  C(N-1,1)= dt/2/h^3;
  C(N,1)= -dt/h^3-dt/h*a(1);
  C(N,2)= dt/2/h^3;
  D = speye(N);
 
  
end
    

    
    



