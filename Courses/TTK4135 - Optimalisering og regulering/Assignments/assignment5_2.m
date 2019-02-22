%% MPC
% System matrices 
A = [   
        0       0       0       ;
        0       0       1       ;
        0.1     -0.79   1.78    
];
b = [1; 0; 0.1];
c = [0 0 1];

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

% Inequality constraints
xlb = -Inf(N*nx,1);    % Lower bound on x
xub =  Inf(N*nx,1);    % Upper bound on x
ulb = -ones(N*nu,1);   % Lower bound on u
uub =  ones(N*nu,1);   % Upper bound on u
lb = [xlb; ulb];      % Lower bound on z
ub = [xub; uub];      % Upper bound on z

X = zeros(nx,N+1);
U = zeros(nu,N);
X(:,1) = [0 0 1]'; % Initial state
for k=1:N

    beq = [A*X(:,k);zeros(nx*(N-1),1)];

    % Solving optimization problem using quadprog
    [z,fval,exitflag,output,lambda] = quadprog(G,[],[],[],Aeq,beq,lb,ub);
    x = [x0(3);z(nx:nx:N*nx)];
    u = z(N*nx+1:N*(nx+nu));
    
    U(k) = u(1);
    
    A = [
            0       0       0;e
            0       0       1;
            0.3     -0.955  1.95
    ];
    b = [1;0;0];
    
    X(:,k+1) = A*X(:,k)+b*U(:,k);
    
end
% Plot optimal trajectory and optimal input
hold on;
subplot(2,1,1);
plot(X(3,:)); % Plot on 0 to N
grid('on');
ylabel('x_t')
hold on;
subplot(2,1,2);
plot(U(1,:)); % Plot on 0 to N
grid('on');
ylabel('u_t');