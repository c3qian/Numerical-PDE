    v = 16;
    x0 = 0;    
    dx = 0.1;
    h = dx;
    N = 16 / dx;
    u_max = -8;
    dt = dx / ((4 / dx ^ 2) + 6 * abs(u_max)) * 0.5;

    t_end = 2;
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
    

    % initialization
    for i = 1 : n_methods
       % part b one soliton
       v_old(:, i) = - 8 ./ cosh(2 * x) .^ 2;
       
       % part c two soliton solition
       %v_old(:,i) = -8./(cosh(2*x +8).^2) -2./(cosh(x).^2);
       
       % part d gaussian profile
       %v_old(:,i) = -5 * exp(-(x./2).^2);
    end
    
    
    % the main iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the first iteration    
    for j = 3 : N - 2
        temp(j, 1) = v_old(j, 1) + dt * ((v_old(j - 1, 1) + ...
            + v_old(j, 1) + v_old(j + 1, 1)) * (v_old(j + 1, 1) - v_old(j - 1, 1)) / (dx) + ...
            -(v_old(j + 2, 1) - 2 * v_old(j + 1, 1) + 2 * v_old(j - 1, 1) - v_old(j - 2, 1)) / (2 * dx ^ 3) );
    end
    
    temp(1, 1) = v_old(1, 1) + dt * ((v_old(N, 1) + v_old(1, 1) + ...
            + v_old(2, 1)) * (v_old(2, 1) - v_old(N, 1)) / (dx) -(v_old(3, 1) - 2 * v_old(2, 1) +...
            + 2 * v_old(N, 1) - v_old(N - 1, 1)) / (2 * dx ^ 3) );
    temp(2, 1) = v_old(2, 1) + dt * ((v_old(1, 1) + v_old(2, 1) + ...
            + v_old(3, 1)) * (v_old(3, 1) - v_old(1, 1)) / (dx) -(v_old(4, 1) +...
            - 2 * v_old(3, 1) + 2 * v_old(1, 1) - v_old(N, 1)) / (2 * dx ^ 3) );
    temp(N - 1, 1) = v_old(N - 1, 1) + dt * ((v_old(N - 2, 1) + ...
            +v_old(N - 1, 1) + v_old(N, 1)) * (v_old(N, 1) - v_old(N - 2, 1)) / (dx) + ...
            -(v_old(1, 1) - 2 * v_old(N, 1) + 2 * v_old(N - 2, 1) - v_old(N - 3, 1)) / (2 * dx ^ 3) );
    temp(N, 1) = v_old(N, 1) + dt * ((v_old(N - 1, 1) + v_old(N, 1) +...
            + v_old(1, 1)) * (v_old(1, 1) - v_old(N - 1, 1)) / (dx) -(v_old(2, 1) +...
            - 2 * v_old(1, 1) + 2 * v_old(N - 1, 1) - v_old(N - 2, 1)) / (2 * dx ^ 3) );
       
    % after the first iteration
    for iter = 2 : n_it 
        % method 1: leapfrog scheme
        for j = 3 : N - 2
            v_new(j, 1) = v_old(j, 1) + 2 * dt * ((temp(j - 1, 1) + temp(j, 1) + temp(j + 1, 1)) * (temp(j + 1, 1) - temp(j - 1, 1)) / (dx) - (temp(j + 2, 1) - 2 * temp(j + 1, 1) + 2 * temp(j - 1, 1) - temp(j - 2, 1)) / (2 * dx ^ 3) );
        end
        
        v_new(1, 1) = v_old(1, 1) + 2 * dt * ((temp(N, 1) + temp(1, 1) + temp(2, 1)) * (temp(2, 1) - temp(N, 1)) / (dx - (temp(3, 1) - 2 * temp(2, 1) + 2 * temp(N, 1) - temp(N - 1, 1)) / (2 * dx ^ 3) ));
        v_new(2, 1) = v_old(2, 1) + 2 * dt * ((temp(1, 1) + temp(2, 1) + temp(3, 1)) * (temp(3, 1) - temp(1, 1)) / (dx)+...
            - (temp(4, 1) - 2 * temp(3, 1) + 2 * temp(1, 1) - temp(N, 1)) / (2 * dx ^ 3) );
        v_new(N - 1, 1) = v_old(N - 1, 1) + 2 * dt * ((temp(N - 2, 1) + temp(N - 1, 1) + temp(N, 1)) * (temp(N, 1) + ...
            - temp(N - 2, 1)) / (dx) - (temp(1, 1) - 2 * temp(N, 1) + 2 * temp(N - 2, 1) - temp(N - 3, 1)) / (2 * dx ^ 3) );
        v_new(N, 1) = v_old(N, 1) + 2 * dt * ((temp(N - 1, 1) + temp(N, 1) + temp(1, 1)) * (temp(1, 1) - temp(N - 1, 1)) / (dx)+ ...
            - (temp(2, 1) - 2 * temp(1, 1) + 2 * temp(N - 1, 1) - temp(N - 2, 1)) / (2 * dx ^ 3) );
        
        % method 2: Goda scheme
        % since we update the diagnal(-2)(2) at every iteration, so we
        % simply put it in the main loop.
        e = ones(N-1,1);
        e1 = ones(N,1);
        for i = 2:N
            dia1(i-1,1) = temp(i,2) + temp(i-1,2);     
        end

        for i = 1: N-1
            dia1(i,2) = temp(i,2) + temp(i+1,2);
        end
        %dia1(N,1) = 0;
        %dia1(N,2) = 0;
        % some parameter to initialize the matrix C.
      para1 = -0.5*dt/h^3*e1;
      para2 = [dt/h^3 * e + dt/h * dia1(1:N-1,1);0];
      para3 = [0; -dt/h^3 - dt/h * dia1(1:N-1,2)];
      para4 = 0.5*dt/h^3*e1; 

      C = spdiags([para1 para2 e1 para3 para4], [-2 -1 0 1 2], N, N);
        % upper right and lower left
      C(1,N) = dt/h^3 - dt/h * dia1(1,1);
      C(1,N-1) = -0.5*dt/h^3;
      C(2,N) = -0.5*dt/h^3;
      C(N,1) = -dt/h^3 + dt/h * dia1(N,2);
      C(N-1,1) = 0.5*dt/h^3;
      C(N,2) = 0.5*dt/h^3;
      D = speye(N);
        
        v_new(:,2)=C\(D*v_old(:,2));
        temp(:,2) = v_new(:,2);
        % the exact solution given by part one.      
        %v_exact(:)=-8./cosh(2*(x-16*iter*dt)).^2;
        v_exact(:) = -0.5*v./(cosh(0.5*sqrt(v)*(x-v*t(iter+1)-x0))).^2;
        
        % graphical output      
        plot(x, v_exact(:),'^r-')         
        hold on       
        for i=1:n_methods
            plot(x, v_new(:,1), '*g-')
            plot(x, v_new(:,2), '*b-')
        end
        legend('exact','LF','Goda','FontSize',6);
        hold off       
        %pause(0.001)
        v_old = temp;
        temp = v_new;
        
    end
  
