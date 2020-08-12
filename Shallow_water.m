%Assignment 3

function Shallow_water
    g=1;
    N=100;
    h=2/N;
    dimension=3;
    x=(-1+h/2:h:1-h/2);
    %v_new=zeros(N,dimension);
    v_old=zeros(N,dimension);
    
    t_end=0.5;
    t = 0;

    % initial conditions
    v_old(:,1)=2 .*(x<0 & x>-1)+(x>=0 & x<=1).*1;
    v_old(:,2)=0 .*(x<0 & x>-1)+(x>=0 & x<=1).*0;
    v_old(:,3)=0 .*(x<0 & x>-1)+(x>=0 & x<=1).*0;
 
    while t <=t_end
        % first column is h
        % second column is hu
        % third column is hv
        dt=0.5* min(h./(abs(v_old(:,2)./v_old(:,1))+sqrt(g.*v_old(:,1))));    
        % method 1: LF
        % using ghost cells as following
        aug1 = [[v_old(2:N,1);1] [v_old(2:N,2);0] [v_old(2:N,3);0]];
        aug2 = [[2;v_old(1:N-1,1)] [0;v_old(1:N-1,2)] [0;v_old(1:N-1,3)]];
        v_new =v_old -dt/h*(LF(v_old ,aug1)-LF(aug2,v_old));   
        v_old= v_new;
        t = t+dt;
    end
  
    hold on
    
    plot(x,v_new(:,3)./v_new(:,1),'+g-')
    plot(x,v_new(:,2)./v_new(:,1),'*R-')
    plot(x,v_new(:,1),'*b-')
    
    axis([-1 1 -1 3])
    xlabel('x')
    ylabel('v')
    title('Shallow Water simulation t=0.5')
    legend('LF','location','southeast')
    
    hold off
end
 
% LF flux function
function F=LF(v_left,v_right)
    g=1;
    lambda_left = abs(v_left(:,2)./v_left(:,1))+sqrt(g.*v_left(:,1));
    lambda_right = abs(v_right(:,2)./v_right(:,1))+sqrt(g.*v_right(:,1));
    lambda_max=max(lambda_left,lambda_right)/2;
    F(:,1)=(v_left(:,2)+v_right(:,2))/2-lambda_max.*(v_right(:,1)-v_left(:,1));
    F(:,2)=(v_left(:,2).^2./v_left(:,1)+g.*v_left(:,1).^2./2 + v_right(:,2).^2./v_right(:,1)+g.*v_right(:,1).^2./2)/2-lambda_max.*(v_right(:,2)-v_left(:,2));
    F(:,3)=(v_left(:,3).*v_left(:,2)./v_left(:,1)+v_right(:,3).*v_right(:,2)./v_right(:,1))/2-lambda_max.*(v_right(:,3)-v_left(:,3));
end
 

 
