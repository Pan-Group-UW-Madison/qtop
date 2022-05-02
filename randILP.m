clear all;
close all;
clc;

N = 5000;
M = 5;

n = 0;

for i = 1:M
    f = 100 * rand(1, N) - 50;
    lb = zeros(1, N);
    ub = ones(1, N);
    intcon = 1:N;
%     Aeq = randi(N/2, 1, N);
    Aeq = ones(1, N);
    beq = N/2.3;
    x1 = intlinprog(f, intcon, [], [], Aeq, beq, lb, ub);
    x2 = linprog(f, [], [], Aeq, beq, lb, ub);

    if norm(x1 - x2, 1) > 0
        n = n + 1;
    end

end

disp(n);
