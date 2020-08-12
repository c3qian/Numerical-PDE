%Assignment 3

function Dam_break
    g=1;
    N=60;
    hx=2/N;
    hy=2/N;
    dimension=3;
    x=(-1+hx/2:hx:1-hx/2);
    y=(-1+hy/2:hy:1-hy/2);
    [X,Y] = meshgrid(x,y);
    %v_new=zeros(N,dimension);
    v_old=zeros(N,N,dimension);
    
    t_end=3;
    t = 0.1;
%     % problem 1
%     a=1;
%     b=2;
%  initial condition

    %v_old(:,:,1)= 2*(x<=0.5 & x>=-0.5)*ones(length
%     v_old(:,:,2)=0.*(x<=0.5 & x>=-0.5)'*(y<=0.5 & y>=-0.5)+(x<=1 & x>=-1)'*(y>0.5 & y<1).*0 ...
%         + (x<=1 & x>=-1)'*(y>-1 & y<-0.5).*0 +(x<-0.5& x>-1)'*(y>-0.5&y<0.5).*0 +(x<1& x>0.5)'*2 (y>-0.5&y<0.5).*0;
%     v_old(:,:,3)=0.*(x<=0.5 & x>=-0.5)'*(y<=0.5 & y>=-0.5)+(x<=1 & x>=-1)'*(y>0.5 & y<1).*0 ...
%         + (x<=1 & x>=-1)'*(y>-1 & y<-0.5).*0 +(x<-0.5& x>-1)'*(y>-0.5&y<0.5).*0 +(x<1& x>0.5)'* (y>-0.5&y<0.5).*0;
    v_old(:,:,2) = zeros(N);
    v_old(:,:,3) = zeros(N);
    while t <= t_end
        dx_lamx=min(hx./(abs(v_old(:,:,2)./v_old(:,:,1))+sqrt(g.*v_old(:,:,1))),[],'all');
        dy_lamy=min(hy./(abs(v_old(:,:,3)./v_old(:,:,1))+sqrt(g.*v_old(:,:,1))),[],'all');
        dxy_max = min(dx_lamx,dy_lamy);
        dt= 0.4 * dxy_max;   
        % method 1: LF
        aug1 = cat(3,[v_old(2:N,:,1);ones(N,1)'],[v_old(2:N,:,2);zeros(N,1)'],[v_old(2:N,:,3);zeros(N,1)']);
        aug2 = cat(3,[2.*ones(N,1)';v_old(1:N-1,:,1)],[zeros(N,1)';v_old(1:N-1,:,2)],[zeros(N,1)';v_old(1:N-1,:,3)]);
        aug3 = cat(3,[v_old(:,2:N,1),-1*ones(N,1)],[v_old(:,2:N,2),zeros(N,1)],[v_old(:,2:N,3),zeros(N,1)]);
        aug4 = cat(3,[-2*ones(N,1),v_old(:,1:N-1,1)],[zeros(N,1),v_old(:,1:N-1,2)],[zeros(N,1),v_old(:,1:N-1,3)]);
        v_new = v_old -(dt/hx).*(LF(v_old,aug1)-LF(aug2,v_old))-(dt/hy).*(LFG(v_old,aug3)-LFG(aug4,v_old));   
        
        t = t + dt
        v_old= v_new;
    end
    % Problem 1
    % v_exact=(x<t_end)+2*(x>1+2*t_end)+(x>=t_end & x<=1+2*t_end).*((x-t_end)/(1+t_end)+1);
   
    hold on
    colormap winter
    %plot(x,v_exact,'^r-')
    mesh(X,Y,real(v_new(:,:,1)))
    %plot(x,v_new(:,1),'*b-')
    %plot(x,v_new(:,2),'+g-')
    %plot(x,v_new(:,3),'+g-')
    %axis([-1 1 -1 1 0 3])
    xlabel('x')
    ylabel('y')
    zlabel('v')
    title('Shallow Water simulation')
    %legend('LF','location','southeast')
    
    hold off
end
 
% % LF flux function
% function F =LF(v_left,v_right)
%     g=1;
%     lambda_left = abs(v_left(:,2)./v_left(:,1))+sqrt(g.*v_left(:,1));
%     lambda_right = abs(v_right(:,2)./v_right(:,1))+sqrt(g.*v_right(:,1));
%     lambda_max=max(lambda_left,lambda_right)/2;
%     F(:,1)=(v_left(:,2)+v_right(:,2))/2-lambda_max.*(v_right(:,1)-v_left(:,1));
%     F(:,2)=(v_left(:,2).^2./v_left(:,1)+g.*v_left(:,1).^2./2 ...
%         + v_right(:,2).^2./v_right(:,1)+g.*v_right(:,1).^2./2)/2-lambda_max.*(v_right(:,2)-v_left(:,2));
%     F(:,3)=(v_left(:,3).*v_left(:,2)./v_left(:,1)...
%         +v_right(:,3).*v_right(:,2)./v_right(:,1))/2-lambda_max.*(v_right(:,3)-v_left(:,3));    
% end
% function G = LFG(v_left,v_right)
%     g=1;
%     lambda_left = abs(v_left(:,3)./v_left(:,1))+sqrt(g.*v_left(:,1));
%     lambda_right = abs(v_right(:,3)./v_right(:,1))+sqrt(g.*v_right(:,1));
%     lambda_max=max(lambda_left,lambda_right)/2;
%     G(:,1)=(v_left(:,2)+v_right(:,2))/2-lambda_max.*(v_right(:,1)-v_left(:,1));
%     G(:,3)=(v_left(:,2).^2./v_left(:,1)+g.*v_left(:,1).^2./2 ...
%         + v_right(:,2).^2./v_right(:,1)+g.*v_right(:,1).^2./2)/2-lambda_max.*(v_right(:,2)-v_left(:,2));
%     G(:,2)=(v_left(:,3).*v_left(:,2)./v_left(:,1)...
%         +v_right(:,3).*v_right(:,2)./v_right(:,1))/2-lambda_max.*(v_right(:,3)-v_left(:,3));
% end
% LF flux function
function F =LF(v_left,v_right)
    g=1;
    lambda_left = min(abs(v_left(:,:,2)./v_left(:,:,1))+sqrt(g.*v_left(:,:,1)),[],'all');
    lambda_right = min(abs(v_right(:,:,2)./v_right(:,:,1))+sqrt(g.*v_right(:,:,1)),[],'all');
    lambda_max=max(lambda_left,lambda_right)/2;
    F(:,:,1)=(v_left(:,:,2)+v_right(:,:,2))/2-lambda_max.*(v_right(:,:,1)-v_left(:,:,1));
    F(:,:,2)=(v_left(:,:,2).^2./v_left(:,:,1)+g.*v_left(:,:,1).^2./2 ...
        + v_right(:,:,2).^2./v_right(:,:,1)+g.*v_right(:,:,1).^2./2)/2-lambda_max.*(v_right(:,:,2)-v_left(:,:,2));
    F(:,:,3)=(v_left(:,:,3).*v_left(:,:,2)./v_left(:,:,1)...
        +v_right(:,:,3).*v_right(:,:,2)./v_right(:,:,1))/2-lambda_max.*(v_right(:,:,3)-v_left(:,:,3));    
end
function G = LFG(v_up,v_down)
    g=1;
    lambda_left = abs(v_up(:,:,3)./v_up(:,:,1))+sqrt(g.*v_up(:,:,1));
    lambda_right = abs(v_down(:,:,3)./v_down(:,:,1))+sqrt(g.*v_down(:,:,1));
    lambda_max=max(lambda_left,lambda_right)./2;
    G(:,:,1)=(v_up(:,:,2)+v_down(:,:,2))/2-lambda_max.*(v_down(:,:,1)-v_up(:,:,1));
    G(:,:,3)=(v_up(:,:,2).^2./v_up(:,:,1)+g.*v_up(:,:,1).^2./2 ...
        + v_down(:,:,2).^2./v_down(:,:,1)+g.*v_down(:,:,1).^2./2)/2-lambda_max.*(v_down(:,:,2)-v_up(:,:,2));
    G(:,:,2)=(v_up(:,:,3).*v_up(:,:,2)./v_up(:,:,1)...
        +v_down(:,:,3).*v_down(:,:,2)./v_down(:,:,1))/2-lambda_max.*(v_down(:,:,3)-v_up(:,:,3));
end
