    MEM_LAT = 100;
    L3_LAT = 10;
    L2_LAT = 5;
    L1_LAT = 1;
    L = 4;
    N = 1024;
    m_c = N;
    n_c = N;
    k_c = N;

    % Define the objective function (Memory Access Time)
    objective = @(x, m,n,k) ( ...
          (m*k)*MEM_LAT ...                  
        + ((m*n*k)/x(1))*MEM_LAT ...
        + ((m*n*k)/(x(3)))*L3_LAT ...
        + (m*k*n*(1/x(4)))*L2_LAT ...
        + 2*(m*k*n*(1/x(2)))*MEM_LAT ...
        + (m*k*n*(1/x(5)))*L1_LAT ...
        + (m*n*k*(1/x(4)))*L2_LAT );

    % Initialize sum
%     total_sum = (2*N*N*N)/objective([x1 x2 x3 x4 x5], N,N,N);
    total_sum = 0;
    syms x1 x2 x3 x4 x5
    START = 512;
    END = N;
    % Loop over m from 512 to 1026
    for m = START:END
        % Calculate objective function for current m
%         We are minimizing total time not total q, because reducing
%         latency is the ultimate measure
%         total_sum = total_sum + (2*m*m*m)/objective([x1 x2 x3 x4 x5], m,m,m);
        total_sum = total_sum + objective([x1 x2 x3 x4 x5], m,m,m);
    end
%total_sum
% total_sum = -1*total_sum; % allow maximization so fmincon would minimize total_sum, and maximize -1*totalsum
total_sum_f = matlabFunction(total_sum);
obj_f = @(x) total_sum_f(x(1), x(2), x(3), x(4), x(5));

occupancy = 0.9;
% Define the nonlinear constraint function
nonlcon = @(x) deal([
     x(1) * x(2) - 4194304*occupancy ;   % a*b < 4194304
     x(3) * x(2) - 131072*occupancy;    % c*b < 131072
     x(4) * x(2) + x(2)*x(5) - 8192*occupancy;      % d*b < 8192
     x(4) * x(5) - 64;         % d*e < 32
%      x(1) - m_c;
%      x(2) - k_c;
%      x(3) - n_c;

], []); % Empty array for equality constraints
% Define the lower bounds
lb = [0, 0, 0, 16, 4];

% Define the upper bounds (assuming no upper bounds)
ub = [N, N , N, 16, 4 ];

% Provide an initial guess
x0 = [1, 1, 1, 16, 4];

% Solve the optimization problem
[x, fval] = fmincon(obj_f, x0, [], [], [], [], lb, ub, nonlcon);

% Display the optimal solution
disp('Optimal solution:');
disp(['MC = ', num2str(x(1))]);
disp(['KC = ', num2str(x(2))]);
disp(['NC = ', num2str(x(3))]);
disp(['MR = ', num2str(x(4))]);
disp(['NR = ', num2str(x(5))]);
disp(['Objective value = ', num2str(fval)]);

