% Winter 2020
% Assignment 2

% first initialize some parameters  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 70;
% spatial step
h =7/N;
umax = 2;
dt = 0.5 * h/umax;
n_methods = 2;
t_end = 2;
n_it = t_end/dt;
% initialize some arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatial grid 
x=(-1 +h/2:h: 6 -h/2);
% temporal grid
t=(0:dt:n_it*dt);
% arrays for the numerical approximations at times n and n+1
v_new=zeros(N,n_methods);
v_old=zeros(N,n_methods);
alpha=zeros(N,2);
f_star=zeros(N,2);
v_exact= zeros(N,n_methods);


for i = 1:n_methods
    % initial profile
    %problem 1:
    a = 1;
    b = 2;
    v_old(:,i) = 1*(x < 0) + (x>=0 & x<=1).*(1 + x) + 2*(x > 1);
    
% %     % Problem 2:
%      a = 2;
%      b = 1;
%      v_old(:,i)=2 * (x<0) + (x >= 0 & x <= 1).*(2 - x) + 1*(x > 1);
%     % Bonus quesion
%      a = 0;
%      b = 0;
%      v_old(:,i)= 0 * (x<0) + (x >= 0 & x <= 1).* x + 0 *(x > 1);
end

for iteration = 1:n_it
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % method 1: 
    % problem 1:   
    % calculate f(j+1/2)
    % append two ghost cells at the beginning and the back
     alpha(1:N,1) = max(abs(v_old(1:N,1)),abs([v_old(2:N,1);b]));
%     % calculate f(j-1/2)
     alpha(1:N,2) = max(abs([a;v_old(1:N-1,1)]),abs(v_old(1:N,1)));
     f_star(1:N,1) = (v_old(1:N,1).^2/2 + [v_old(2:N,1);b].^2/2)/2 - alpha(1:N,1)./2 .*([v_old(2:N,1);b]-v_old(1:N,1));
     f_star(1:N,2) = ([a;v_old(1:N-1,1)].^2./2 + v_old(1:N,1).^2./2)/2 - alpha(1:N,2)./2 .*(v_old(1:N,1) - [a;v_old(1:N-1,1)]);
     v_new(1:N,1) = v_old(1:N,1)-dt/h.*(f_star(1:N,1) - f_star(1:N,2));
    
%     %problem 2:
%     %everything is the same except the value of a and b
%      alpha(1:N,1) = max(abs(v_old(1:N,1)),abs([v_old(2:N,1);b]));
%     %calculate f(j-1/2)
%      alpha(1:N,2) = max(abs([a;v_old(1:N-1,1)]),abs(v_old(1:N,1)));
%      f_star(1:N,1) = (v_old(1:N,1).^2/2 + [v_old(2:N,1);b].^2/2)/2 - alpha(1:N,1)./2 .*([v_old(2:N,1);b]-v_old(1:N,1));
%      f_star(1:N,2) = ([a;v_old(1:N-1,1)].^2./2 + v_old(1:N,1).^2./2)/2 - alpha(1:N,2)./2 .*(v_old(1:N,1) - [a;v_old(1:N-1,1)]);
%      v_new(1:N,1) = v_old(1:N,1)-dt/h.*(f_star(1:N,1) - f_star(1:N,2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % method 2:
    % first step
    [v_plus_half_left,v_plus_half_right,v_minus_half_left,v_minus_half_right] = SLR(v_old(:,2),N,a,b);
    v_n_half = v_old(:,2)-dt/2 .*(NF(v_plus_half_left,v_plus_half_right)- NF(v_minus_half_left,v_minus_half_right))/h;
    % second step
    [v_plus_half_left,v_plus_half_right,v_minus_half_left,v_minus_half_right] = SLR(v_n_half,N,a,b);
    v_new(:,2) = v_old(:,2)-dt*(NF(v_plus_half_left,v_plus_half_right)- NF(v_minus_half_left,v_minus_half_right))/h;
    v_old = v_new;
end
% the exact solution
% problem 1
 v_exact(:,1) = (x < t_end)+(x >= t_end & x <= 1 + 2*t_end).*((x-t_end)/(1 + t_end)+1) + 2*(x > 1 + 2 * t_end);

% Problem 2
% Because there contains a shock, so we have to divide into two cases based
% on t_end

 if t_end < 1
     v_exact(:,2) = 2 * (x < 2*t_end) + (x > 1 + t_end)+ +( x >= 2 * t_end & x <= 1+t_end).*((x-2 * t_end)/(t_end - 1) + 2);
 else
     v_exact(:,2) = 2 * (x < (t_end-1) * 1.5 + 2)+(x > (t_end - 1)* 1.5 + 2)+ 1.5*( x==(t_end - 1) * 1.5 + 2);
 end
 
 % Bonous question
 v_exact(:,3) = 0 * (x < 0) + x/(1+t_end) .* (x>= 0 & x <= sqrt(1 + t_end)) + 0* (x > sqrt(1 + t_end));


hold on
plot(x,v_exact(:,1),'^r-')
plot(x,v_new(:,1),'*g-')
plot(x,v_new(:,2),'+b-')
axis([1 6 0.5 2.5])
xlabel('x')
ylabel('v')
title('1D Burgers Equation P1')
%title('Bonus question')
legend('exact','Lax-F','second accuracy','location','southeast')
%legend('exact','second accuracy','location','southeast')
%legend('exact','Lax-F','location','southeast')
hold off


% Numerical flux function
function f = NF(v_minus,v_plus)
    Beta = max(abs(v_minus),abs(v_plus));
    f = (v_minus.^2./2 + v_plus.^2./2)/2 - Beta./2.* (v_plus - v_minus);
end
 
%second order using linear reconstruction
% V_i-1 append a to the left
% V_i+1 append b to the right
% continue to add ghost cells for the case V_i-2, V_i+2
% VL represents the Value of linear construction
function [v_plus_half_left,v_plus_half_right,v_minus_half_left,v_minus_half_right] = SLR(VL,N,a,b)
    r = (VL-[a;VL(1:N-1)])./([VL(2:N);b]- VL);
    r_left = ([a;VL(1:N-1)]-[a;a;VL(1:N-2)])./(VL- [a;VL(1:N-1)]);
    r_right = ([VL(2:N);b]-VL)./([VL(3:N);b;b]-[VL(2:N);b]);
    % v_i+1/2,left
    v_plus_half_left = VL+ 1/2* max(0,min(r,1)).*([VL(2:N);b]- VL);
    v_plus_half_right = [VL(2:N);b]- 1/2* max(0,min(r_right,1)).*([VL(3:N);b;b]- [VL(2:N);b]);
    % v_i-1/2,left
    v_minus_half_left = [a;VL(1:N-1)]+ 1/2* max(0,min(r_left,1)).*(VL-[a;VL(1:N-1)]);
    v_minus_half_right = VL - 1/2*max(0,min(r,1)).*([VL(2:N);b]- VL);
end


