% MIT License

% Copyright (c) 2018 Anguigems

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Differences over time for 2 Double Pendulum's simultaneously.
% ODE solver using RK4 method. 
% t=time, y=theta
% compares difference in angle theta w.r.t time

domain = [0,100];       % Domain to solve ODE of order k
Narray=[1000];         % All the number of steps to try. If this is an array, the code will plot them all
y0=[0.4, 0, 2, 0, 0.417, 0, 2, 0]; % Initial Conditions:[theta1 velocity1 theta2 velocity2] [y1 y1' y2 y2']
% both pendulums display different solutions for even 1 degree=0.017
% radians different.
[m1,m2,m3,m4,l1,l2,l3,l4,g,u1,u2] = param(); % calling parameters
% the initial value of y up to its (k-1)th derivatives- possibly an array
% Output: t, solution  y


for i =1:length(Narray)
    N=Narray(i);
    [t,y]= eulerODE(domain,y0,N);
    figure('Name','Differences Over Time','NumberTitle','off')
    plot(t,abs(abs(y(:,1))-abs(y(:,5))),'-r', 'Linewidth',2.5, 'DisplayName','upper masses') % upper mass m1.     
    hold on
    plot(t,abs(abs(y(:,3))-abs(y(:,7))),'-b','Linewidth',2.5, 'DisplayName','lower masses') % lower mass m2.
    hold on
    legend('upper masses','lower masses', 'Location', 'northwest') 
end

% Cosmetics for graphs
set(gca,'FontSize',14)
ylim([0,1.5])
xlim([0,40])
xlabel('time t (seconds)')
ylabel('angle \theta (radians)')
title('Difference Over Time')


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
% The current step is determined by the previous step added to the average of the previous step and the estimation of the current step.  
end


function f = dydt(t,y)
% Put the RHS of your ODE (system) here.
[m1,m2,m3,m4,l1,l2,l3,l4,g,u1,u2] = param();
% f and y are arrays
f(1) = y(2);
f(2) = (-g*sin(y(1))-u1*l2*y(4)^2*sin(y(1)-y(3))+u1*g*sin(y(3))*cos(y(1)-y(3))-u1*l1*y(2)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l1*(1-u1*(cos(y(1)-y(3)))^2));
f(3) = y(4);
f(4) = (-g*sin(y(3))+l1*y(2)^2*sin(y(1)-y(3))+g*sin(y(1))*cos(y(1)-y(3))+u1*l2*y(4)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l2*(1-u1*(cos(y(1)-y(3)))^2));
% The first f is column 2 of y? Yes it is.
f(5) = y(6);
f(6) = (-g*sin(y(5))-u2*l4*y(8)^2*sin(y(5)-y(7))+u2*g*sin(y(7))*cos(y(5)-y(7))-u2*l3*y(6)^2*sin(y(5)-y(7))*cos(y(5)-y(7)))/(l3*(1-u2*(cos(y(5)-y(7)))^2));
f(7) = y(8);
f(8) = (-g*sin(y(7))+l3*y(6)^2*sin(y(5)-y(7))+g*sin(y(5))*cos(y(5)-y(7))+u2*l4*y(8)^2*sin(y(5)-y(7))*cos(y(5)-y(7)))/(l4*(1-u2*(cos(y(5)-y(7)))^2));
end


% Put the physical parameters here.
function [m1,m2,m3,m4,l1,l2,l3,l4,g,u1,u2] = param()
m1=1; m2=1; m3=1; m4=1; l1=0.5; l2=0.5; l3=0.5; l4=0.5; g=9.8; u1=m2/(m1+m2); u2=m4/(m3+m4);
% l1 is the length of the rod between the pivot and bob 1.
% l2 is the length of the rod between bob 1 and bob 2.
% g is gravity.
end
