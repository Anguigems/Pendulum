% ODE solver using RK4 method for a Double Pendulum.
% t=time, y=theta

domain = [0,30];       % Domain to solve ODE of order k
Narray=[1000];         % All the number of steps to try. If this is an array, the code will plot them all
y0=[0.2 0 2 0]; % Initial Conditions:[theta1 velocity1 theta2 velocity2]
[m1,m2,l1,l2,g,u] = param(); % calling parameters
% the initial value of y up to its (k-1)th derivatives- possibly an array
% Output: t, solution  y


for i =1:length(Narray)
    N=Narray(i);
    [t,y]= eulerODE(domain,y0,N);
    figure('Name','Double Pendulum RK4 Solution','NumberTitle','off')
    plot(t,y(:,1),'--b', 'Linewidth', 1.5, 'DisplayName','upper mass') % upper mass m1
    hold on
    plot(t,y(:,3),'--r','Linewidth',1.5, 'DisplayName','lower mass') %lower mass m2
    hold on
%    Try changing y(:,1) to y(:,2) to plot velocity rather than angle theta.
%    and changing y(:,3) to y(:,4) to plot velocity rather than angle theta.
%    and change y limits to at least [-20,20]. 
    legend('upper mass','lower mass')
end


% Cosmetics for graphs
set(gca,'FontSize',16)
%legend('Trapezium','Simpson', 'Location', 'northwest')
ylim([-3,3])
xlim([0,15])
xlabel('time t (seconds)')
ylabel('angle \theta (radians)')
title('Double Pendulum RK4 Solution')


function [t,y]= eulerODE(domain,y0,N)
% Euler's Method for Solving IVP dydt = f(t,y)
% y and f are possibly column vectors (i.e. arrays)
% Input: domain, initial value y0 (an array), number of steps N
% Output:  t and y for each step
t(1)=domain(1); 
y(1,:)=y0; % tells it to set row 1 as y0 vector.
h=(domain(2)-domain(1))/N;         % this is the step size
for i=1:N
  t(i+1)=t(i)+h;
  y(i+1,:)=ODEstep(t(i),y(i,:),h);
% above line states i+1th row = ODEstep of ith row. 
end
end


function y=ODEstep(t,y,h)
% RK4 method
s_1=dydt(t,y);
s_2=dydt(t+(h/2),y+(h/2)*s_1);
s_3=dydt(t+(h/2),y+(h/2)*s_2);
s_4=dydt(t+h,y+h*s_3);
y=y+(h/6)*(s_1+2*(s_2+s_3)+s_4);
% The current step is determined by the previous step added to the 
% average of the previous step and the estimation of the current step.  
end


function f = dydt(t,y)
% Put the RHS of your ODE (system) here.
[m1,m2,l1,l2,g,u] = param();
% f and y are arrays
f(1) = y(2);
f(2) = (-g*sin(y(1))-u*l2*y(4)^2*sin(y(1)-y(3))+u*g*sin(y(3))*cos(y(1)-y(3))-u*l1*y(2)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l1*(1-u*(cos(y(1)-y(3)))^2));
f(3) = y(4);
f(4) = (-g*sin(y(3))+l1*y(2)^2*sin(y(1)-y(3))+g*sin(y(1))*cos(y(1)-y(3))+u*l2*y(4)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l2*(1-u*(cos(y(1)-y(3)))^2));
end


% Put the physical parameters here.
function [m1,m2,l1,l2,g,u] = param()
m1=1; m2=1; l1=0.6; l2=0.3; g= 9.8; u= m2/(m1+m2);
% l1 is the length of the rod between the pivot and bob 1.
% l2 is the length of the rod between bob 1 and bob 2.
% m1 is the upper mass.
% m2 is the lower mass.
% g is gravity.
end

