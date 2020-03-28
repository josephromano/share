close all
x=linspace(0,1,100);
A=0.5;
phi=pi/2;
y1 = A*sin(2*pi*x);
y2 = A*sin(2*pi*x+phi);
plot(x,y1,'k',x,y2,'k','linewidth',2)
grid on
xlim([-0.1 1.1])
ylim([-1 1])
xlabel('x (m)','fontsize',14)
ylabel('y (m)','fontsize',14)
%axis equal
print -djpeg superposition.jpg
