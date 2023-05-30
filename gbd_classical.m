clear all;
close all;
clc;

nelx = 120;
nely = 40;

display = 0;

figure;
tiledlayout(1,1, 'TileSpacing', 'none', 'Padding', 'none')
set(gcf, 'position', [200, 200, 200*ceil(nelx/nely), 200])

optTime = 0;
t1 = tic;
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
KE = 1 / (1 - nu ^ 2) / 24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11]);
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
%% PREPARE FILTER
iH = ones(nelx * nely * (2 * (ceil(rmin) - 1) + 1) ^ 2, 1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;

for i1 = 1:nelx

    for j1 = 1:nely
        e1 = (i1 - 1) * nely + j1;

        for i2 = max(i1 - (ceil(rmin) - 1), 1):min(i1 + (ceil(rmin) - 1), nelx)

            for j2 = max(j1 - (ceil(rmin) - 1), 1):min(j1 + (ceil(rmin) - 1), nely)
                e2 = (i2 - 1) * nely + j2;
                k = k + 1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0, rmin - sqrt((i1 - i2) ^ 2 + (j1 - j2) ^ 2));
            end

        end

    end

end

H = sparse(iH, jH, sH);
Hs = sum(H, 2);
%% INITIALIZE ITERATION
x = ones(nely, nelx);
loop = 0;

dvol = 200;
stage = 1;
totalLoop = 0;
epsilon = 5e-4;

numFEM = 0;
numFeasibleCut = 0;

compliance = [];
volume = [];

while stage < 3
    loop = loop + 1;
    Lower = 0;
    Upper = 1e9;

    if stage == 1
        vol = floor(volfrac * nelx * nely);
        vol = vol - dvol;
        vol = max(vol, volfrac0 * nelx * nely);
        volfrac = vol / (nelx * nely);

        if (abs(vol - volfrac0 * nelx * nely) < 1)
            epsilon = 5e-4;
        end

    end

    if stage == 2
        vol = floor(volfrac0 * nelx * nely);
        epsilon = 5e-4;
        break;
    end

    if display == 1
        set(gca, 'Position', [0, 0, 1, 1])
        colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;
    end

    innerLoop = 0;
    exitFlag = 1;
    xOptimal = x;
    xTarget = [];
    ceTarget = [];
    cTarget = [];

    xFeasible = [];
    ceFeasible = [];
    cFeasible = [];

    % prepare the intitial solution
    if (loop == 1)
        xPhys = x;
        xPhys(xPhys < 0.01) = Emin;
        sK = reshape(KE(:) * (Emin + xPhys(:)' * (E0 - Emin)), 64 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
        c = sum(sum((Emin + xPhys * (E0 - Emin)) .* ce));
        ce(:) = ce .* xPhys;
        ce(:) = H * (ce(:) ./ Hs);

        numFEM = numFEM + 1;
    end

    optTimer1 = tic;
    [xResult, cost, exitFlag] = gbdMasterCut(reshape(x, [], 1), c, reshape(ce, 1, []), [], [], [], vol);
    optTime = optTime + toc(optTimer1);

    x = reshape(xResult, size(x, 1), size(x, 2));
    xOptimal = x;

    if display == 1
        set(gca, 'Position', [0, 0, 1, 1])
        colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;
    end

    while (1)
        innerLoop = innerLoop + 1;

        % primal problem
        xPhys = x;
        xPhys(xPhys < 0.01) = Emin;
        sK = reshape(KE(:) * (Emin + xPhys(:)' * (E0 - Emin)), 64 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
        c = sum(sum((Emin + xPhys * (E0 - Emin)) .* ce));
        ce(:) = ce .* xPhys;
        ce(:) = H * (ce(:) ./ Hs);

        numFEM = numFEM + 1;

        if c < Upper
            xOptimal = x;
            Upper = c;
        end

        xTarget = [xTarget reshape(x, [], 1)];
        ceTarget = [ceTarget; reshape(ce, 1, [])];
        cTarget = [cTarget; c];
        index = [];

        for i = 1:length(cTarget)

            if (cTarget(i) <= c)
                index = [index; i];
            end

        end

        % master problem
        optTimer1 = tic;
        [xResult, cost, exitFlag] = gbdMasterCut(xTarget(:, index), cTarget(index), ceTarget(index, :), xFeasible, cFeasible, ceFeasible, vol);
        optTime = optTime + toc(optTimer1);

        if exitFlag == 1
            x = reshape(xResult, size(x, 1), size(x, 2));
        else
            break;
        end

        if display == 1
            set(gca, 'Position', [0, 0, 1, 1])
            colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;
        end

        if cost > Upper || abs(Upper - cost) / Upper < epsilon
            compliance = [compliance; Upper];
            volume = [volume; vol];

            fprintf(' It.:%5i index: %5i Obj.:%11.4f Vol.:%7.3f Gap.:%7.3f%%\n', loop, length(index), Upper, mean(x(:)), (Upper - cost) / Upper * 100);

            break;
        end

        compliance = [compliance; Upper];
        volume = [volume; vol];

        fprintf(' It.:%5i index: %5i Obj.:%11.4f Vol.:%7.3f Gap.:%7.3f%%\n', loop, length(index), Upper, mean(x(:)), (Upper - cost) / Upper * 100);

    end

    x = xOptimal;

    totalLoop = totalLoop + innerLoop;

    if stage == 1 && volfrac <= volfrac0
        stage = 2;
    elseif stage == 2
        stage = 3;
    end

end

totalTime = toc(t1);
disp(['total time: ', num2str(totalTime), 's']);
disp(['optimization time: ', num2str(optTime), 's']);
disp(['num of fem: ', num2str(numFEM)]);

set(gca, 'Position', [0, 0, 1, 1])
colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

function [x, objFunc, exitFlag] = gbdMasterCut(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    n = size(y, 2);

    if n > 1
        l = size(y, 1);
        f = zeros(1, l + 1);
        f(1) = 1;
        A = zeros(n + 1, l + 1);
        b = zeros(n + 1, 1);

        for i = 1:n
            A(i, 1) = -1;
            A(i, 2:end) = -weight(i, :);
            b(i) = -obj(i) -weight(i, :) * y(:, i);
        end

        A(n + 1, :) = ones(1, l + 1);
        A(n + 1, 1) = 0;
        b(n + 1) = vol;

        model.obj = f;
        model.A = sparse(A);
        model.rhs = b;

        for i = 1:n
            model.sense(i) = '<';
        end

        model.sense(n + 1) = '=';

        for i = 1:l
            model.vtype(i + 1) = 'B';
        end

        model.vtype(1) = 'C';
        model.lb = zeros(1, l + 1);
        model.lb(1) = -inf;
        model.ub = ones(1, l + 1);
        model.ub(1) = inf;
        model.modelsense = 'min';

        params.outputflag = 1;

        result = gurobi(model, params);

        x = result.x(2:end);
        objFunc = result.x(1);
        exitFlag = 1;

    else
        l = size(y, 1);
        model.obj = -weight(1, :);
        model.A = sparse(ones(1, l));
        model.rhs = vol;
        model.sense = '=';
        model.vtype = 'B';
        model.modelsense = 'min';

        params.outputflag = 0;

        result = gurobi(model, params);

        x = result.x;
        objFunc = obj(1) - weight(1, :) * (x - y(:, 1));
        exitFlag = 1;
    end

end
