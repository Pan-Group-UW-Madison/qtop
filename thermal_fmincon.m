close all;
clear all;
clc;

global nelx nely k0 kmin
nelx = 40;
nely = 40;
k0 = 1;
kmin = 1e-3;

volfrac0 = 0.5;
volfrac = 1.0;

options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true);
x = fmincon(@fem, 0.5*ones(nelx * nely, 1), ones(1, nelx * nely), volfrac0 * nelx * nely, [], [], zeros(1, nelx * nely), ones(1, nelx * nely), [], options);
x = reshape(x, nely, nelx);
colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

function [c, ce] = fem(xPhys)
    global nelx nely k0 kmin
    KE = [2/3 -1/6 -1/3 -1/6
        -1/6 2/3 -1/6 -1/3
        -1/3 -1/6 2/3 -1/6
        -1/6 -1/3 -1/6 2/3];
    nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
    edofVec = reshape(nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
    edofMat = repmat(edofVec, 1, 4) + repmat([0 nely + [1 0] -1], nelx * nely, 1);
    iK = reshape(kron(edofMat, ones(4, 1))', 16 * nelx * nely, 1);
    jK = reshape(kron(edofMat, ones(1, 4))', 16 * nelx * nely, 1);
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F = sparse(1:(nelx + 1) * (nely + 1), 1, 0.01, (nelx + 1) * (nely + 1), 1);
    % F = sparse(ceil((nelx+1)*(nely+1)/2), 1, 1, (nelx + 1) * (nely + 1), 1);
    U = zeros((nely + 1) * (nelx + 1), 1);
    fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
    % fixeddofs = 1:nely+1;
    alldofs = 1:(nely + 1) * (nelx + 1);
    freedofs = setdiff(alldofs, fixeddofs);

    xPhys = reshape(xPhys, nelx, nely);
    sK = reshape(KE(:) * (kmin + xPhys(:)' * (k0 - kmin)), 16 * nelx * nely, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
    c = sum(sum((kmin + xPhys * (k0 - kmin)) .* ce));
    ce = -ce(:);
end
