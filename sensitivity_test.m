clear all;
close all;
clc;

nelx = 30;
nely = 10;
rmin = 2;
volfrac0 = 0.5;
volfrac = 1.0;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = 1 / (1 - nu^2) / 24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11]);
nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(2 * nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] -2 -1], nelx * nely, 1);
iK = reshape(kron(edofMat, ones(8, 1))', 64 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * nelx * nely, 1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2, 1, -1, 2 * (nely + 1) * (nelx + 1), 1);
U = zeros(2 * (nely + 1) * (nelx + 1), 1);
fixeddofs = union([1:2:2 * (nely + 1)], [2 * (nelx + 1) * (nely + 1)]);
alldofs = [1:2 * (nely + 1) * (nelx + 1)];
freedofs = setdiff(alldofs, fixeddofs);

x = ones(nely, nelx);

sK = reshape(KE(:) * (Emin + x(:)' * (E0 - Emin)), 64 * nelx * nely, 1);
K = sparse(iK, jK, sK); K = (K + K') / 2;
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
ce_base = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
obj = F' * U;

ce = zeros(nely, nelx);

for i = 1:nely

    for j = 1:nelx
        x = ones(nely, nelx);
        x(i, j) = 0;
        sK = reshape(KE(:) * (Emin + x(:)' * (E0 - Emin)), 64 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        obj2 = F' * U;
        ce(i, j) = obj2 - obj;
    end

end

difference = ce - ce_base;
