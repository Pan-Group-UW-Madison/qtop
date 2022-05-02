clear all;
close all;
clc;

figure(1);
set(gcf, 'position', [0, 0, 400, 400]);

hold on;
plot([0, 0], [0.45, 0.55], 'r-', 'linewidth', 2);
plot([1, 1], [0.45, 0.55], 'g-', 'linewidth', 2);

plot([0, 0], [0, 0.45], 'k-', 'linewidth', 2);
plot([0, 0], [0.55, 1], 'k-', 'linewidth', 2);
plot([1, 1], [0, 0.45], 'k-', 'linewidth', 2);
plot([1, 1], [0.55, 1], 'k-', 'linewidth', 2);
plot([0, 1], [1, 1], 'k-', 'linewidth', 2);
plot([0, 1], [0, 0], 'k-', 'linewidth', 2);

plot([0.5, 0.5], [-0.05, 1.05], 'b-.', 'linewidth', 2);
plot([-0.05, 1.05], [0.5, 0.5], 'y--', 'linewidth', 2);

xlabel('$x$', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', 20);

axis equal;
