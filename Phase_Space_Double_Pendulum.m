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


% ODE solver using RK4 method to approximate the solution to the Double 
% Pendulum in order to plot Phase Space (Poncaire) diagram.

% t=time, y=theta.

domain = [0,11.98];       % Domain to solve ODE of order k
Narray=[1000];         % All the number of steps to try. If this is an array, the code will plot them all
y0=[pi/2 0 pi/2+0.06 0]; % Initial Conditions:[theta1 velocity1 theta2 velocity2]
[m1,m2,l1,l2,g,u] = param(); % calling parameters
% the initial value of y up to its (k-1)th derivatives- possibly an array
% Output: t, solution  y

% The following loop plots the phase space for each respective mass/bob on
% the double pendulum.
for i =1:length(Narray)
    N=Narray(i);
    [t,y]= eulerODE(domain,y0,N);
    figure('Name','Double Pendulum Phase Space','NumberTitle','off')
    x = mod(y+pi,2*pi); % adjusted values to prevent it going past pi,-pi.
    plot(x(:,1)-pi,y(:,2),'-b', 'Linewidth', 1.5, 'DisplayName','upper mass') % upper mass m1 (angle,velocity) plot
    hold on
    plot(x(:,3)-pi,y(:,4),'-r','Linewidth',1.5, 'DisplayName','lower mass') %lower mass m2 (angle,velocity) plot
end


% Cosmetics for graphs produced by the above.
set(gca,'FontSize',16)
legend('upper mass','lower mass','Location', 'northwest')
ylim([-30,30])
xlim([-3.2,3.2])
xlabel('Angle \theta (radians)')
ylabel('Velocity (radians)')
title('Double Pendulum Phase Space')


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
% RK4 method for solving the ode
s_1=dydt(t,y);
s_2=dydt(t+(h/2),y+(h/2)*s_1);
s_3=dydt(t+(h/2),y+(h/2)*s_2);
s_4=dydt(t+h,y+h*s_3);
y=y+(h/6)*(s_1+2*(s_2+s_3)+s_4);
end


function f = dydt(t,y)
% The ode for the double pendulum has been converted into a linear system
% of equations that can be solved using the RK4 method to gain a numerical
% approximation to the solution for theta (y in the code).
[m1,m2,l1,l2,g,u] = param();
% f and y are arrays
f(1) = y(2);
f(2) = (-g*sin(y(1))-u*l2*y(4)^2*sin(y(1)-y(3))+u*g*sin(y(3))*cos(y(1)-y(3))-u*l1*y(2)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l1*(1-u*(cos(y(1)-y(3)))^2));
f(3) = y(4);
f(4) = (-g*sin(y(3))+l1*y(2)^2*sin(y(1)-y(3))+g*sin(y(1))*cos(y(1)-y(3))+u*l2*y(4)^2*sin(y(1)-y(3))*cos(y(1)-y(3)))/(l2*(1-u*(cos(y(1)-y(3)))^2));
end


% Physical parameters are here, these paraeters can be varied to affect the
% resuting motion of the double pendulum.
function [m1,m2,l1,l2,g,u] = param()
m1=0.09702; 
m2=0.08119; 
l1=0.254; 
l2=0.2286; 
g= 9.8; 
u= m2/(m1+m2);
% l1 is the length of the rod between the pivot and bob 1.
% l2 is the length of the rod between bob 1 and bob 2.
% m1 is the upper mass (the bob between the two beams/rods).
% m2 is the lower mass (the bob at the end of the second beam/rod).
% g is gravity.
end

