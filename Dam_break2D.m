function Dam_break2D
    g=1;
    N=60;
    hx=2/N;
    hy=2/N;
    %dimension=3;
    x=(-1+hx/2:hx:1-hx/2);
    y=(-1+hy/2:hy:1-hy/2);
    [X,Y] = meshgrid(x,y);
    %v_new=zeros(N,dimension);
    v_old_1=zeros(N,N);
    v_old_2=zeros(N,N);
    v_old_3=zeros(N,N);
    t_end=3;
    t = 0;

    
    % initial conditions
    
    for i=1:N
        for j=1:N
            if (x(i)>=-0.5 && x(i)<=0.5) && (y(j)>=-0.5 && y(j)<=0.5)
                    %h
                    v_old_1(i,j)=2;                  
                    %hu
                    v_old_2(i,j)=0;
                    %hv
                    v_old_3(i,j)=0;
            else
                    %h
                    v_old_1(i,j)=1;                  
                    %hu
                    v_old_2(i,j)=0;
                    %hv
                    v_old_3(i,j)=0;             
            end
        end
    end
    
    while t <= t_end
        dx_lamx=min(hx./(abs(v_old_2(:,:)./v_old_1(:,:))+sqrt(g.*v_old_1(:,:))),[],'all');
        dy_lamy=min(hy./(abs(v_old_3(:,:)./v_old_1(:,:))+sqrt(g.*v_old_1(:,:))),[],'all');
        dxy_max = min(dx_lamx,dy_lamy);
        dt= 0.4 * dxy_max;   
        % method 1: LF
        % adding ghost cells
        aug11 = [v_old_1(2:N,:);v_old_1(N,:)];
        aug12 = [v_old_2(2:N,:);-1.*v_old_2(N,:)];
        aug13 = [v_old_3(2:N,:);1.*v_old_3(N,:)];
        
        aug21 = [v_old_1(1,:);v_old_1(1:N-1,:)];
        aug22 = [-1*v_old_2(1,:);v_old_2(1:N-1,:)];
        aug23 = [1*v_old_3(1,:);v_old_3(1:N-1,:)];
        
        augy11 = [v_old_1(:,2:N),v_old_1(:,N)];
        augy12 = [v_old_2(:,2:N),1.*v_old_2(:,N)];
        augy13 = [v_old_3(:,2:N),-1.*v_old_3(:,N)];
        
        augy21 = [v_old_1(:,1),v_old_1(:,1:N-1)];
        augy22 = [1*v_old_2(:,1),v_old_2(:,1:N-1)];
        augy23 = [-1*v_old_3(:,1),v_old_3(:,1:N-1)];
        
        [FL1,FL2,FL3]=LFF(v_old_1,v_old_2,v_old_3,aug11,aug12,aug13);
        [FR1,FR2,FR3]=LFF(aug21,aug22,aug23,v_old_1,v_old_2,v_old_3);
        [GL1,GL2,GL3]=LFG(v_old_1,v_old_2,v_old_3,augy11,augy12,augy13);
        [GR1,GR2,GR3]=LFG(augy21,augy22,augy23,v_old_1,v_old_2,v_old_3);
        
        
        v_new_1 = v_old_1 -(dt/hx).*(FL1-FR1)-(dt/hy).*(GL1-GR1);   
        v_new_2 = v_old_2 -(dt/hx).*(FL2-FR2)-(dt/hy).*(GL2-GR2);
        v_new_3 = v_old_3 -(dt/hx).*(FL3-FR3)-(dt/hy).*(GL3-GR3);
        t = t + dt;
        v_old_1= v_new_1;
        v_old_2= v_new_2;
        v_old_3= v_new_3;
        % check the volume is 5 all the time, i.e., the water is all in the
        % tank
        
        vol = sum(sum(v_new_1(:,:)))*hx*hy

       

        colormap winter      
        mesh(X,Y,v_new_1(:,:))
        axis([-1 1 -1 1 0 3])
        
        xlabel('x')
        ylabel('y')
        zlabel('h')
        title('2D Shallow Water simulation t= 3')
        pause(0.01)
%         hold on
%         plot(t,vol,'+r-')
%         title('total Water volume inside the box')
%         axis([0 3 4 6])
%         hold off
        

    end
end

   
    
    % LF flux function
function [F1,F2,F3] =LFF(v_left_1,v_left_2,v_left_3,v_right_1,v_right_2,v_right_3)
    g = 1;
    lambda_left = min(abs(v_left_2(:,:)./v_left_1(:,:))+sqrt(g.*v_left_1(:,:)),[],'all');
    lambda_right = min(abs(v_right_2(:,:)./v_right_1(:,:))+sqrt(g.*v_right_1(:,:)),[],'all');
    lambda_max=max(lambda_left,lambda_right)/2;
    F1=(v_left_2(:,:)+v_right_2(:,:))./2-lambda_max.*(v_right_1(:,:)-v_left_1(:,:));
    F2=(v_left_2(:,:).^2./v_left_1(:,:)+ g.*v_left_1(:,:).^2./2 ...
        + v_right_2(:,:).^2./v_right_1(:,:)+ g.*v_right_1(:,:).^2./2)./2-lambda_max.*(v_right_2(:,:)-v_left_2(:,:));
    F3=(v_left_3(:,:).*v_left_2(:,:)./v_left_1(:,:)...
        +v_right_3(:,:).*v_right_2(:,:)./v_right_1(:,:))./2-lambda_max.*(v_right_3(:,:)-v_left_3(:,:));    
end

function [G1,G2,G3] = LFG(v_left_1,v_left_2,v_left_3,v_right_1,v_right_2,v_right_3)
    g = 1;
    lambda_left = min(abs(v_left_3(:,:)./v_left_1(:,:))+sqrt(g.*v_left_1(:,:)),[],'all');
    lambda_right = min(abs(v_right_3(:,:)./v_right_1(:,:))+sqrt(g.*v_right_1(:,:)),[],'all');
    lambda_max=max(lambda_left,lambda_right)/2;
    G1=(v_left_3(:,:)+v_right_3(:,:))./2-lambda_max.*(v_right_1(:,:)-v_left_1(:,:));
    G2=(v_left_3(:,:).*v_left_2(:,:)./v_left_1(:,:)...
        +v_right_3(:,:).*v_right_2(:,:)./v_right_1(:,:))./2-lambda_max.*(v_right_2(:,:)-v_left_2(:,:));   
    G3=(v_left_3(:,:).^2./v_left_1(:,:)+ g.*v_left_1(:,:).^2./2 ...
        + v_right_3(:,:).^2./v_right_1(:,:)+ g.*v_right_1(:,:).^2./2)./2-lambda_max.*(v_right_3(:,:)-v_left_3(:,:));
end

