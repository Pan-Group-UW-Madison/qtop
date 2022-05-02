close all;
clear all;
clc;

for i = 1:100
    N = 5000;
    n = 2;
    obj = rand(1, n);
    c = randn(n, N);

    f = zeros(1, N + 1);
    f(1) = 1;
    lb = zeros(1, N + 1);
    lb(1) = -inf;
    ub = ones(1, N + 1);
    ub(1) = inf;

    A = zeros(n + 1, N + 1);
    b = zeros(n + 1, 1);

    for i = 1:n
        A(i, 1) = -1;
        A(i, 2:end) = -c(i, :);
        b(i) = -obj(i);
    end

    intcon = 1:N;
    intcon = intcon + 1;

    A(end, :) = ones(1, N + 1);
    A(end, 1) = 0;
    b(end) = N / 2;

    [x, objFunc, exitFlag, ~] = intlinprog(f, intcon, A, b, [], [], lb, ub);
end
