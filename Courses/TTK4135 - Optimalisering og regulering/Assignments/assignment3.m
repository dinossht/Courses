%% Contour plot with constraints
x1 = linspace(-1,7,7);
x2 = linspace(-1,7,7);
[X,Y] = meshgrid(x1,x2);
Z = -3*X-Y;
contour(X,Y,Z)
hold on;
plot(x1,5-x1/3);
hold on;
plot(x1,8-2*x1);
axis([-1 7 -1 7]);

%% Simplex method to calculate the optimal A and B

clear('all');
clc;

% Problem data: (Remember that all vectors are column vectors.)
c = [-3 -2 0 0]';
A = [2 1 1 0; 1 3 0 1];
b = [8 15]';
% A fesible starting point:
x0 = [0 0 8 15]'; % The starting point must also have four variables!

% Solve problem and print a report: (Call simplex(c,A,b,x0,[]) if you do 
% not want the report printed).
[x, fval, iterates] = simplex(c,A,b,x0,'report');

% Extract iterates in the space of x_1 and x_2:
iter_x1_x2 = iterates(1:2, :);

% Extract iterates as individual vectors (containing x_1 and x_2):
iter_1 = iter_x1_x2(:,1);
iter_2 = iter_x1_x2(:,2);
iter_3 = iter_x1_x2(:,3);


