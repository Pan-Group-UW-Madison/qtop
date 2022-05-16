clc;
kmin = 1e-4;

nelx = 40;
nely = 40;

KE = [2/3 -1/6 -1/3 -1/6
    -1/6 2/3 -1/6 -1/3
    -1/3 -1/6 2/3 -1/6
    -1/6 -1/3 -1/6 2/3];

nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 4) + repmat([0 nely + [1 0] -1], nelx * nely, 1);
iK = reshape(kron(edofMat, ones(4, 1))', 16 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 4))', 16 * nelx * nely, 1);

F = sparse(1:(nelx + 1) * (nely + 1), 1, 1 / nelx / nely, (nelx + 1) * (nely + 1), 1);
U = zeros((nely + 1) * (nelx + 1), 1);

fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
alldofs = 1:(nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);

x = ones(nelx, nelx);
% x(:, 30:40) = 0;
x(36:40, :) = 0;
x(1:5, :) = 0;
x(:, 31:40) = 0;
x(40, 40) = 0;

xPhys = x;
xPhys(xPhys < 1e-3) = kmin;
sK = reshape(KE(:) * (xPhys(:)'), 16 * nelx * nely, 1);
K = sparse(iK, jK, sK); K = (K + K') / 2;
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
ce1 = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
c = F'*U;

disp(ce1(36, 1));
disp(ce1(36, 1)*xPhys(36, 1)^2);
disp(c);
disp('------------');

nelx = 80;
nely = 80;

nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 4) + repmat([0 nely + [1 0] -1], nelx * nely, 1);
iK = reshape(kron(edofMat, ones(4, 1))', 16 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 4))', 16 * nelx * nely, 1);

F = sparse(1:(nelx + 1) * (nely + 1), 1, 1 / nelx / nely, (nelx + 1) * (nely + 1), 1);
U = zeros((nely + 1) * (nelx + 1), 1);

fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
alldofs = 1:(nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);

x = ones(nelx, nelx);
% x(:, 30:40) = 0;
x(71:80, :) = 0;
x(1:10, :) = 0;
x(:, 61:80) = 0;
x(80, 80) = 0;

xPhys = x;
xPhys(xPhys < 1e-3) = kmin;
sK = reshape(KE(:) * (xPhys(:)'), 16 * nelx * nely, 1);
K = sparse(iK, jK, sK); K = (K + K') / 2;
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
ce2 = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
c = F'*U;

disp(ce2(71, 1));
disp(ce2(71, 1)*xPhys(71, 1)^2);
disp(c);
disp('------------');

nelx = 160;
nely = 160;

nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 4) + repmat([0 nely + [1 0] -1], nelx * nely, 1);
iK = reshape(kron(edofMat, ones(4, 1))', 16 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 4))', 16 * nelx * nely, 1);

F = sparse(1:(nelx + 1) * (nely + 1), 1, 1 / nelx / nely, (nelx + 1) * (nely + 1), 1);
U = zeros((nely + 1) * (nelx + 1), 1);

fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
alldofs = 1:(nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);

x = ones(nelx, nelx);
% x(:, 30:40) = 0;
x(141:160, :) = 0;
x(1:20, :) = 0;
x(:, 121:160) = 0;
x(160, 160) = 0;

xPhys = x;
xPhys(xPhys < 1e-3) = kmin;
sK = reshape(KE(:) * (xPhys(:)'), 16 * nelx * nely, 1);
K = sparse(iK, jK, sK); K = (K + K') / 2;
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
ce3 = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
c = F'*U;

disp(ce3(141, 1));
disp(ce3(141, 1)*xPhys(141, 1)^2);
disp(c);