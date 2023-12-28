% Numerical methods course in physics LAB

% Created by Shoaib Mustafa
% SP21-BPH-068
% 6TH SEMESTER PHYSICS COMSATS
% You are free to modify or copy.
% generated for my own practice, sharing so others can benefit. 
% you are encouraged to read the accompanying pdf's 
% each section can be run by 'ctrl/command + enter'


%lab 20 practice 
%ode45


clc;
close all; 
clear all; 



%% Lagrange Interpolation 
data_x = [1, -4, 0]; 

data_y = [3, 13, -23]; 

basis_poly(1) = 0;

syms x; 


for i = 1:length(data_x)

    x_i = data_x(i);

    my_prod = 1; 

    %loops over all data points. 
    for j =1:length(data_x)

        %for each x_i it will loop over all the data
        %and add to product. except when i = j

        x_j = data_x(j);

        if j ~= i

            my_prod = my_prod * (x - x_j) / (x_i - x_j)
        
        end

    end
    basis_poly(i) = my_prod;
end

lagrange_poly = 0; 


for k=1:length(data_x)
    lagrange_poly = lagrange_poly + basis_poly(k) * data_y(k);
end

%simplifies the polynomial. adds subtracts etc. 
lagrange_poly = simplify(lagrange_poly);

fprintf('Lagrange polynomial is %s', lagrange_poly)


%% Interp1 easier interpolation 

%this function returns an interpolated function. easily. 
clear all
clear variables
close all

data_x = [1, -4, 0]; 

data_y = [3, 13, -23]; 

interp1(data_x, data_y, 4)
%interp2 gives you data for two variable system. 
%interp3 gives you three variable. 

% Problem: interpolate surface. 
[X, Y] = meshgrid(-2:2)
R = sqrt(X.^2 + Y.^2);
V = sin(R) ./ R; 

[Xq, Yq] = meshgrid(-2: 0.2 : 2);

ans2 = interp2(X, Y, V, Xq, Yq);

surf(X, Y, V)
figure
surf(Xq, Yq, ans2)



%% Interpolation with Spline 

 
clear all; 
clear variables; 
close all; 

data_x = [-3, -2, -1, 0, 1,2,3];
data_y = [-1, -1, -1, 0, 1,1,1]; 

xq = -3:0.1:3; 

ans = spline(data_x, data_y, xq)

plot(data_x, data_y, 'b')
hold on 
plot(ans)
axis square
hold offnote



%% Numerical Integration Trapz
% if y = [c1 , c2 , c2 ] trapz(y) => [trapz(c1), trapz(c2), trapz(c3)]
%Q = trapz(Y) %unit spacing

% use trapz to calculate integral I = lim([1 5]) (x^2 dx) 
clear all;
clear variables; 

x_interval = 1: 5; 
%x_interval = linspace(1, 5); 
f = x_interval .^ 2;

trapz_ans = trapz(f); 

% double integration. NEST 

% syntax
% trapz(interval, function, column/row)


% Quad 
% uses adaptive simpson quadrature
% quad(function, lim_a, lim_b, tol)
% dblquad(fun, xmin, xmax, ymin, ymax) for double integration. 
% fun should be matrix values. 
% cannot handle infinity. 

% integral 
% integral(fun, xmin, xmax) 
% can handle infinity. 


% integral2 and integral 3
% integral2 has same syntax. 
% integral3 also same but three things 
% integrand has to be function handle. 



% int function 
% computes indefinite integral with symbolic variable. 
% no inf 
% syntax int(function, variable to integrate, lim_a, lim_b);
syms x;
f = sin(x); 
int_int = int(f); 
% can calculate pole of a function by 
% int(f, x, 0, 2, 'PrincipalValue', true);




%% 1D radioactive decay
% by Kevin Berwick,
% based on 'Computational Physics' book by N Giordano and H Nakanishi
% Section 1.2 p2
% Solve the Equation dN/dt = -N/tau
N_uranium_initial = 1000; %initial number of uranium atoms
npoints = 100; %Discretize time into 100 intervals
dt = 1e7; % time step in years
tau=4.4e9; % mean lifetime of 238 U
N_uranium = zeros(npoints,1); % initializes N_uranium, a vector of dimension npoints X 1,to being
all zeros
time = zeros(npoints,1); % this initializes the vector time to being all zeros
N_uranium(1) = N_uranium_initial; 
% the initial condition, first entry in the vector N_uranium is N_uranium_initial
time(1) = 0; %Initialise time
for step=1:npoints-1 % loop over the timesteps and calculate the numerical solution
N_uranium(step+1) = N_uranium(step) - (N_uranium(step)/tau)*dt;
time(step+1) = time(step) + dt;
end
% For comparison , calculate analytical solution
t=0:1e8:10e9;
N_analytical=N_uranium_initial*exp(-t/tau);
% Plot both numerical and analytical solution
plot(time,N_uranium,'r')
hold on
plot(t,N_analytical,'b'); %plots the numerical solution in red and the analytical solution in blue
xlabel('Time in years');
ylabel('Number of atoms');




%% ODE eular 
%given that dy/dt + 20y = 7*exp(-0.5*t) compute y(0.2) by taking h= 0.1  
% y(0) = 5
% know that tn = t0 + n*h 
% y0 = 5 
% t0 = 0 
% update y0 and t0 in each iter
%define step size h = b-a / n
% f = input('enter your functoin: '); 

f = @(t,y) -20*y + 7*exp(0.5*t);
t0 = 0; 
y0 = 5; 
h = 0.02;
tn = 0.2 % point at which you'd like to evaluate and stop 
%for iterations 
% n = b - a / h 

n = (tn - t0)/h; 

%matlab can't evaluate t at zero 
%so giving initial conditions for the for loop 

t(1) = t0; 
y(1) = y0; 

for i =1:n
    y(i+1) = y(i) + h*f(t(i), y(i));
    t(i+1) = t0 + i*h; 
    fprintf('y(%.2f) = %.4f\n', t(i+1), y(i+1))
end




%% RK method 2
%general equation of runge kutta. 
% k1 = hf(tn, yn)
% k2 = hf(tn + h, yn+k1);
% and yn+1 = yn + 1/2[k1 + k2] 

%given that dy/dy + 20*y = 7*exp(-0.5*t);
%compute y(0.2) using RK method of order 2. by taking h = 0.1; 
clc;
clear all;
close all; 

f = @(t,y) -20*y + 7*exp(0.5*t);
t0 = 0; 
y0 = 5; 
h = 0.1; 
tn = 0.2; 
n = (tn - t0)/ h ; 


y(1) = y0; 
t(1) = t0; 

for i=1:n
    k1 = h*f(t(i),y(i));
    t(i+1) = t0 + i*h; 
    k2 = h*f(t(i+1), y(i) + k1);
    y(i+1) = y(i) + (1/2)*(k1 + k2); 
    fprintf('y(%.2f) = %.4f\n', t(i+1),y(i+1))
end


%% What about rk2 for couple DE?
%lets say we have y''(t) - 0.05y'(t) + 0.15*y(t) = 0 
%funda is simple. if we have a 2nd order DE 
%make substituion to convert to single order de y' = z and y'' = z'
%change the initial conditions accordingly. y'(0) = z(0) and y''(0) = z'(0)
%define two function that are interlinked. 
% dy/dt = z and dz/dt = 0.05z - 0.15y
% note that now our functions will be of f(t,y,z) since interlinked
%define the start and stop points. 
%define the step size. 
%run the loops for n number of steps 
%only addition for
clc; 
clear all;
close all; 

f = @(t,y,z) z;
g = @(t,y,z) 0.05*z - 0.15*y;

t0 = 0; 
y0 = 1; %first dependent variable.
z0 = 0; %2nd dependent variable. 

h = 0.5;
tn = 1; % point at which you'd like to evaluate and stop 
%for iterations 
% n = b - a / h 

n = (tn - t0)/h; 

%matlab can't evaluate t at zero 
%so giving initial conditions for the for loop 

t(1) = t0; 
y(1) = y0;
z(1) = z0; 


%for two equations. 
for i=1:n
    t(i+1) = t0 + i*h; 

    k1_f = h*f(t(i),y(i), z(i));

    k1_g = h*g(t(i),y(i), z(i));

    k2_f = h*f(t(i+1), y(i) + k1_f, z(i) + k1_g);

    %notice the k1 not g
    k2_g = h*g(t(i+1), y(i) + k1_f, z(i) + k1_g);
    
    y(i+1) = y(i) + (1/2)*(k1_f + k2_f); 
    
    z(i+1) = z(i) + 1/2*(k1_g + k2_g); 

    fprintf('y(%.2f) = %.4f\n', t(i+1),y(i+1))
    fprintf('z(%.2f) = %.4f\n', t(i+1),z(i+1))
end

%% RK method 4 
%equations for rk4 is 
% k1 = hf(tn,yn);
% k2 = hf(tn+h/2, yn+ k1/2) 
% k3 = hf(tn+h/2, yn+ k2/2) 
% k4 = hf(tn+h, yn + k3)
% y_n+1 = yn + 1/6 [k1 + 2k2 + 2k3 + k4) ] 

%so similarly. 
clc; 
clear all variables; 
close all; 

f = @(t,y) -20*y + 7*exp(-0.5*t);
t0 = 0;
tn = 0.2; 
y0 = 5; 
h = 0.01; 
n = (tn - t0) / h; 


y(1) = y0;
t(1) = t0; 

for i =1:n
    t(i+1) = t0 + i*h;
    k1 = h*f(t(i), y(i));
    k2 = h*f(t(i) + h/2, y(i) + k1/2);
    k3 = h*f(t(i) + h/2, y(i) + k2/2); 
    k4 = h*f(t(i) + h, y(i) + k3); 
    y(i+1) = y(i) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
    fprintf('y(%.2f) = %.4f\n', t(i+1), y(i+1))
end


%% simple ode
x_i = 0; 
x_f = 10; 


x =[x_i, x_f]
IC=2; 
[X, Y] = ode45(@(x,y) func_1(x,y), x, IC); 
plot(X, Y)
xticks(0:1:10)
yticks(1:.2:3)
ylabel('y', Interpreter='latex')


%% A bit complex. 
IC = 10^-3; 
x_i = 0; 
x_f = 1.5; 

x= [x_i, x_f];

[X, Y] = ode45(@(x,y) func_2(x,y), x, IC)
plot(X,Y,'Marker','o','Color','r'); 
xticks(0:.5:1.5)
yticks(-7:1:2)


%% another one 
x_i = 0; 
x_f = 1; 
x = [x_i, x_f];
IC = 1; 
[X, Y] = ode45(@(x,y) func_3(x,y), x, IC)
plot(X,Y)


%% Some notes 
% Tic toc 
% just a note, tic toc is a timimg function 
% when you call tic it starts the timer 
% when you call toc it stops the timer and prints the time elapsed. 

% diff can do numerical as well as numerical differentiation. 
% diff(f, var, n) is the syntax. 

% gradient(F) returns gradient of one dimension partial f by partial x
% spacing assumed is 1 
% gradient(F, h) 



%% coupled ode
%Through Eular.
%lets say we have y''(t) - 0.05y'(t) + 0.15*y(t) = 0 
%funda is simple. if we have a 2nd order DE 
%make substituion to convert to single order de y' = z and y'' = z'
%change the initial conditions accordingly. y'(0) = z(0) and y''(0) = z'(0)
%define two function that are interlinked. 
% dy/dt = z and dz/dt = 0.05z - 0.15y
% note that now our functions will be of f(t,y,z) since interlinked
%define the start and stop points. 
%define the step size. 
%run the loops for n number of steps 



f = @(t,y,z) z;
g = @(t,y,z) 0.05*z - 0.15*y;

t0 = 0; 
y0 = 1; %first dependent variable.
z0 = 0; %2nd dependent variable. 

h = 0.5;
tn = 1; % point at which you'd like to evaluate and stop 
%for iterations 
% n = b - a / h 

n = (tn - t0)/h; 

%matlab can't evaluate t at zero 
%so giving initial conditions for the for loop 

t(1) = t0; 
y(1) = y0;
z(1) = z0; 

for i =1:n
    %change here from eular is that function havea addition of z in it. 
    y(i+1) = y(i) + h*f(t(i), y(i), z(i));
    z(i+1) = z(i) + h*g(t(i), y(i), z(i)); 
    t(i+1) = t0 + i*h; 
    fprintf('y(%.2f) = %.4f\n', t(i+1), y(i+1))
    fprintf('z(%.2f) = %.4f\n', t(i+1), z(i+1))

end



%% Composite Trapezoid Rule 



%% points inside circle x^2 + y^2 =2  
total_points = 500; 
x = -2 + 4 * rand(total_points,1); 
y = -2 + 4 * rand(total_points,1); 

r = sqrt(x.^2 + y.^2);

pici = find(r<2);
x_inside = x(pici); 
y_inside = y(pici);

plot(x, r)



%% Lab 25 
% Monte carlo methods. 
% A square is generated of length 0<= (x,y) <=1 length r 
% area of square is known. r * r = r^2 
% area of pi is known. pi*r^2
% divide that by 4 to get only one quadrant of pi. 
% pi * r ^2 / 4 
% ratio of circle within square is (pi * r^2 /4) / r^2 ) 
% r^2 gets cancelled. ratio of area is ' pi / 4 '
% that is the area of pi inside one quadrant only. 
% in order to get the result for full pi, we need to multiply the result by
% 4 to include 4 quadrants. 

%first populating the square with random numbers
total_points = 5000;
x = rand(total_points, 1); 
y = rand(total_points, 1); 

%why rand? cause our circle has radius 1. if it had radius 2 we would use
%the a + (b-a) * rand(points , 1) formula. 

% pouring points into the equation to generate set of points that satisfy
% the equation. 

r = sqrt(x.^2 + y.^2); 
%if it was ellipse then? (x.^2 / a^2) + (y.^2 / b^2) = c 
%where a and b represent major and minor axis. considerably harder to plot


%now the step is simple. 
%find points inside of r matrix, which is less than 1. 
%this would mean we are trying to find points that satisfy the inequality
% x^2 + y^2 <= 1 
% any points that satisfy this are deemed to be inside of circle radius 1
% rest of the points are still inside the square but not inside the circle.

points_inside_circle = find( r<= 1);
%this returns the index numbers of such points. not actual values. 
points_outside_circle = find( r> 1); 


% now we have found points inside the circle. 
% but in order to plot we need to find them inside our x and y matrices
% so we simply search for these points in x and y respectively. 

inside_points_x = x(points_inside_circle);
inside_points_y = y(points_inside_circle);

%we are simply telling it to fetch the index numbers values mentioned in that variable
%from x matrix and put that into the left hand side variable. 
%notice the size will not be 5000 anymore. as some will satisfy this and
%some will not. 

%similarly for outside points. 
outside_points_x = x(points_outside_circle); 
outside_points_y = y(points_outside_circle); 

%since monte carlo, we count number of points inside circle. 
%divide that by total number of points. 
%we get ratio that are inside vs that which are inside. 
%that ratio as we said earlier is the area of our chunk of pi.
%multiply that by 4 we get full pi. 
no_points_inside = length(points_inside_circle);

value_of_pi = 4 * no_points_inside/ total_points;
fprintf('Value of pi is %3.4f', value_of_pi);

%also plot 
length(inside_points_x)
length(inside_points_y)
plot(inside_points_x, inside_points_y, 'b'); 
hold on 
plot(outside_points_x, outside_points_y, 'r');
hold off
axis square
%now because x and y were the same lenght. so if (x,y) is a point. that 
%x points will be equal to y points. which helps us in evading an error
%that x and y lengths are not same. scary error that is honestly. 


%what would happend if we were to do the same for ellipse? 
%axis need to be taken care of. lets try? 
%read this to learn a bit about ellipses. https://www.toppr.com/guides/maths/conic-sections/equations-of-ellipse/


%% Generate random points inside ellipse 
% x^2 + 4y^2 = 4 
% rearrange equation. divide by 4. 
% x^2 / 2^2 + y^2 = 1 
% so x axis length -2 to 2 
% y axis -1 to 1. if centered on origin. 

% we simply generate points inside square. then see which points satisfy
% equation. 


total_points = 5000; 
x = -2 + (2 + 2) * rand(total_points, 1); 
y = -1 + (1 + 1) * rand(total_points, 1); 

%using rejection technique. 

r = (x^.2 / 4) + y.^2; 



%% plot Ellipse by MCM
%because area might be a bit complicated. 



%plot can be made more fine if points increased. 
%% lab 26 

%% lab 26 random walk in 1 dimension
% method taught by sir. comparitively easier. 
steps = 100; %num of steps to walk 
probRight = 0.5 %probabillity of moving right in each step


% Initialize position array 
position = zeros(1, steps +1);

for step = 2:steps + 1
    randomNum = rand(); 

    if randomNum < probRight
        position(step) = position(step - 1) + 1; 
    else
        position(step) = position(step -1) - 1; 
    end
end

%plot the random walk 
figure; 
plot(0:steps, position, 'o', 'LineWidth',2);
grid on; 
positiveSteps = position(2:end) > position(1:end-1);
negetiveSteps = position(2:end) < position(1:end-1);

hold on; 
plot(find([0 positiveSteps]), position([1 positiveSteps]), 'go', 'MarkerSize', 8);
plot(find([0 negetiveSteps]), position([1 negetiveSteps]), 'ro', 'MarkerSize', 8);

hold off; 

%% easier random walk 
total_step = 20;
x0 = 0;
y0 = 0;
tx = 0.2;
ty = 0;

plot(x0, y0, 'b')
xlim([-1, 1])
ylim([-0.2, 0.2])
set(gca, 'XTickLabel', [], 'YTickLabel', [])  % Fix: use 'gca' to access the current axes
title('Starting Point')
hold on
pause(1)  % Fix: use parentheses for the pause function

for i = 1:tot_step
    r = rand;
    if r < 0.5
        quiver(x0, y0, tx, ty, 0, 'b')
        
        x0 = x0 + tx;
    else
        quiver(x0, y0, -tx, ty, 0, 'm')
        x0 = x0 - tx;  % Fix: move left, so subtract tx
    end
    hold on
    pause(1)  % Fix: use parentheses for the pause function
end
%%  not taught by professor. harder
A = rand(1,5);
i = 1;
a = 0.5;
ht = .5;
rightt = [];
leftt = [];
while i<=length(A)
    if A(i)<.5
        b = a + .1;
        rightt(i) = A(i);
        annotation('arrow',[b,a],[ht,ht],'color','b')
        text(a,ht,num2str(i))
        a = b;
    else
        b = a - .1;
        leftt(i) = A(i);
        annotation('arrow',[b,a],[ht,ht],'color','r')
        a = b;
        text(a,ht,num2str(i))
    end
    i = i+1;
    pause(2)
end
grid on
disp(A)
disp(rightt)
disp(leftt)


% lab 26 and 27 are almost same. 
% supporting material not present. 

%% lab 27 MCM Integration 
%Question: determine volume of the region
%whose points statisfy inequalities
% 0<= (x,y,z) <= 1 
% x^2 + siny <= z 
% x - z + e^y <= 1 

%concept of mcm integration 
% first generate 'n' random points inside cube of volume 1
% then see how many of them 'm' satisfy last two inequalities
% then the volume of the desired region is approximately m/n 

total_points = 5000; 

%since interval is 0 -> 1, the rand function
% already generate uniform points in this domain
% no extra step is necessary except 

x = rand(total_points, 1); 
y = rand(total_points, 1); 
z = rand(total_points, 1); 

%created x y and z matrices of size 5000 by 1 
% or 5000 rows by 1 column. 

%we just feed values of x, y and z into inequalities
%they will generate a set of points. 

condition_1 = x.^2 + sin(y);
condition_2 = x - z + exp(y); 

%head function helps you peek into first 5 values
%to confirm that you are generating needed values 
%head(condition_1)
%head(condition_2)

%since monte carlo is based on probabillity 
%we count how many occurances of these inequalities 
%is satisfied. we increment the counter. 
%but we have to run a for loop that will check it 5000 times. 

%initialise points counter. 
points_included = 0; 

%do not forget to enclose each condition in a whole bracket 
%and do not add a bracket that encloses whole of if statement. 

for i = 1: total_points
    if (condition_1(i) <= z(i)) && (condition_2(i) <= 1)
        points_included = points_included + 1; 
    end
end 

%for loop will start from 1 till 5000. each increment value will be 
%stored in i variable. 
%we then then i'th index of condition_1 matrix to see if it's less than or
%equal to z. but not any z, i'th z. 
%and condition 2 similarly checking if each index is less than or equal to
%1. if both true. we count. 
%this is to generate number of points inside the inequality. 

reg_vol = points_included / total_points;
fprintf('Approximated volume of the region = %3.4f', reg_vol)


%% Lab 27 Numerical integration using MCM  


%we use formula 
% (multiplication_factor) * summation of all values generated by function
% divided by total number of points. 
% (b-a) * sum(f) / total_points. 
% this one's easy. 
%first define limits as 

x_i = -1;
y_i = -1;
z_i = -1; 
x_f = 1; 
y_f = 1; 
z_f = 1; 

%now create a multiple factor. since we are multiplying all the intervals
%anyway so, (b1_b2)* (a1_a2) etc. 

intervals = (x_f - x_i) * (y_f - y_i) * (z_f - z_i);

%now generate random points in the range provided. 
%remember formula, a + (b-a) * rand(points,1);
%where a is initial and b is final range
%again remember why we aren't using randi in here. because we don't need
%integer values. we need decimal values between interval. 
%and we are giving it explicity range to overcome the limitation of the
%rand function which generates only values between 0 and 1 always. 

x = x_i + (x_f - x_i) * rand(total_points, 1); 
y = y_i + (y_f - y_i) * rand(total_points, 1); 
z = z_i + (z_f - z_i) * rand(total_points, 1); 

%try using head function to see some values generated. 
%now create a matrix which contains calculated values of the integrand. 
f = x.^2 + y.^2 + z.^2; 
%why aren't we defining an anonymous function? as f = @(x,y,z)
%..because we do not need to evaluate the function with new values. 
%we just pump in the ranges of values generated in x y and z and generate a
%single matrix of f which contains 5000 function values

%we still need the summation of the generated values. not the individual
%vlaues or in matrix form. so we use sum 

summed_values = sum(f); 
%takes values from f matrix and sums all of it. returns to that variable. 
%now putting everything into the formula. 

integrated_value = (intervals) * summed_values / total_points; 
%why aren't we doing .*? because we don't have a matrix to evaluate element
%wise or in other words one by one. 

fprintf('Approximated value of the integral= %3.4f', integrated_value); 
%familiarize yourself with fprintf function. its easier. 


%References 
% a good source: https://www.youtube.com/watch?v=_JIupZJEM-k&list=PLO-6jspot8AIxruWnsX6x2g7-Ri2vo64u
% visit links 
% https://en.wikipedia.org/wiki/Monte_Carlo_integration
% watch this video. https://www.youtube.com/watch?v=MKnjsqYVG4Y





