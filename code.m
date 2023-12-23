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



%% coupled ode
%to be filled. 


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

% visit links 
% https://en.wikipedia.org/wiki/Monte_Carlo_integration
% watch this video. https://www.youtube.com/watch?v=MKnjsqYVG4Y





