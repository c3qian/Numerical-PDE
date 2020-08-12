clc;
clear all;
close all;
a = 1;
j = 0;

for j=100:100:1000
 j
 h = 10./j;
 k = 0.5*h;
 n = j+1;
 x = h * [0:j]-5;
 for i=1:n
 uexact(i) = 1 + H(x(i)+1-1) - H(x(i)-1-1);
 end
 for i=1:n
 u(i) = 1 + H(x(i)+1) - H(x(i)-1);
 end
 uo = u;
 for t = 0:k:1
 for i = 2:n-1
 u(i) = 0.5*(uo(i+1)+uo(i-1)) - a*k*(uo(i+1)-uo(i-1))/(2*h);
 end
 uo = u;
 end
 plot(x,u,'--rs','LineWidth',2)
 hold on;
 plot(x,uexact,'-b','LineWidth',2)
 hold off;
 axis([-5 5 0.99 2.01])
 set(gca,'FontSize',30);
 xLabel('x','FontSize',30,'fontweight','b');
 yLabel('u','FontSize',30,'fontweight','b');
 pause(0.01);
 legend('Lax?Friedrichs','Exact');
 error(j/100) = sum(abs(u-uexact))/j;
end
figure;
loglog(10./([100:100:1000]),error,'--r','LineWidth',3);
hold on
set(gca,'FontSize',30);
xLabel('h','FontSize',30,'fontweight','b');
yLabel('error','FontSize',30,'fontweight','b');