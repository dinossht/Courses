%% d)
% System matrices 
A = [   
        0       0       0       ;
        0       0       1       ;
        0.1     -0.79   1.78    
];
b = [1; 0; 0.1];
c = [0 0 1];

x0 = [0 0 1]'; % Initial state

N = 30;  % Length of time horizon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost function 
q = 2*diag([0 0 1]);
r = 2*1;
I_N = eye(N);

R = kron(I_N,r);
Q = kron(I_N,q);

G = blkdiag(Q,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = size(A,2);  % number of states
nu = size(b,2);  % number of inputs
% Equality constaints
Aeqx = zeros(nx*N,nx*N);
Aeqx(nx+1:end,1:end-nx) = kron(eye(N-1),-A);
Aeqx = Aeqx + eye(nx*N);
Aequ = kron(I_N,-b);
Aeq = [Aeqx,Aequ];

beq = [A*x0;zeros(nx*(N-1),1)];

% KKT conditions (eq.16.4)
KKT_mat = [   
        G      -Aeq'                ;
        Aeq     zeros(nx*N)    
];
KKT_vec = [zeros(N*(nx+nu),1);beq];

% KKT solution
KKT_sol = KKT_mat\KKT_vec;
x = [x0(3);KKT_sol(nx:nx:N*nx)];
u2 = KKT_sol(N*nx+1:N*(nx+nu));

% Plot optimal trajectory and optimal input
subplot(2,1,1);
plot(x); % Plot on 0 to N
grid('on');
ylabel('x_t')
subplot(2,1,2);
plot(u2); % Plot on 0 to N
grid('on');
ylabel('u_t');

%% e)
% Solving optimization problem using quadprog
z = quadprog(G,[],[],[],Aeq,beq);
x2 = [x0(3);z(nx:nx:N*nx)];
u2 = z(N*nx+1:N*(nx+nu));

% Plot optimal trajectory and optimal input
subplot(2,1,1);
plot(x2); % Plot on 0 to N
grid('on');
ylabel('x_t')
subplot(2,1,2);
plot(u2); % Plot on 0 to N
grid('on');
ylabel('u_t');

%% f)

% Inequality constraints
xlb = -Inf(N*nx,1);    % Lower bound on x
xub =  Inf(N*nx,1);    % Upper bound on x
ulb = -ones(N*nu,1);   % Lower bound on u
uub =  ones(N*nu,1);   % Upper bound on u
lb = [xlb; ulb];      % Lower bound on z
ub = [xub; uub];      % Upper bound on z

% Solving optimization problem using quadprog
[z,fval,exitflag,output,lambda] = quadprog(G,[],[],[],Aeq,beq,lb,ub);
x3 = [x0(3);z(nx:nx:N*nx)];
u3 = z(N*nx+1:N*(nx+nu));

% Plot optimal trajectory and optimal input
subplot(2,1,1);
plot(x3); % Plot on 0 to N
grid('on');
ylabel('x_t')
subplot(2,1,2);
plot(u3); % Plot on 0 to N
grid('on');
ylabel('u_t');