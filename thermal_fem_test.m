clear all;
close all;
clc;

figure(1);
set(gcf, 'position', [0, 0, 400, 400]);

figure(2);
set(gcf, 'position', [500, 0, 400, 400]);

figure(3);
set(gcf, 'position', [0, 500, 400, 400]);

figure(4);
set(gcf, 'position', [500, 500, 400, 400]);

nel_list = [20, 40, 80, 160, 320, 640];

kmin = 1e-5;
k0 = 1.0;

for i = 1:length(nel_list)
    nelx = nel_list(i);
    nely = nel_list(i);

    h = 1 / nelx;

    KE = 1 / h^2 * [2/3 -1/6 -1/3 -1/6
                -1/6 2/3 -1/6 -1/3
                -1/3 -1/6 2/3 -1/6
                -1/6 -1/3 -1/6 2/3];

    nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
    edofVec = reshape(nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
    edofMat = repmat(edofVec, 1, 4) + repmat([0 nely + [1 0] -1], nelx * nely, 1);
    iK = reshape(kron(edofMat, ones(4, 1))', 16 * nelx * nely, 1);
    jK = reshape(kron(edofMat, ones(1, 4))', 16 * nelx * nely, 1);
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F = ones((nelx + 1) * (nely + 1), 1);
    F(nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20)) = 0;
    F((nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20)) + nelx * (nely + 1)) = 10;
    U = zeros((nely + 1) * (nelx + 1), 1);

    fixeddofs = union(nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20), (nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20)) + nelx * (nely + 1));

    alldofs = 1:(nely + 1) * (nelx + 1);
    freedofs = setdiff(alldofs, fixeddofs);

    x = ones(nely, nelx);
    x(:, ceil(nelx / 2) + 1:nelx) = 1e-3;
    xPhys = x;

    sK = reshape(KE(:) * xPhys(:)', 16 * nelx * nely, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    F(freedofs) = F(freedofs) - K(freedofs, fixeddofs) * F(fixeddofs);
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    U(fixeddofs) = F(fixeddofs);

    figure(1);
    x = linspace(0, 1, nelx + 1);
    y = linspace(0, 1, nely + 1);
    [x, y] = meshgrid(x, y);
    s = surf(x, y, reshape(U, nelx + 1, nely + 1));
    s.EdgeColor = 'none';
    a = colorbar;
    view(0, 90);
    
    xlabel('$x$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$', 'Interpreter', 'latex', 'fontsize', 20);
    
    ylabel(a, '$T$', 'Interpreter', 'latex', 'fontsize', 20);

    figure(2);
    hold on;
    x = linspace(0, 1, nelx + 1);
    T = U(ceil(nely / 2) + 1:nely + 1:(nelx + 1) * (nely +1));
    plot(x, T, 'linewidth', 2);
    
    xlabel('$x$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$T$', 'Interpreter', 'latex', 'fontsize', 20);
    
    figure(3);
    hold on;
%     plot(linspace(0, 1, nelx + 1), U(1:nely+1));
    plot(x, T, 'linewidth', 2);
    xlim([0.45, 0.55]);
    ylim([0 10]);
    
    xlabel('$x$', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('$T$', 'Interpreter', 'latex', 'fontsize', 20);
    
    figure(4);
    hold on;
    plot(linspace(0, 1, nelx + 1), U((nelx+1)*nely+1:end));
end
