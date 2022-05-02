clear all;
close all;
clc;

% figure;
% set(gcf, 'position', [200, 200, 600, 200]);

totalTime = 0;
optTime = 0;

t1 = tic;

nelx = 32;
nely = 3;
nel = nelx * nely + 4;
rmin = 2;
volfrac0 = 98;
volfrac = 100;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = 1 / (1 - nu^2) / 24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11]);
nodenrs = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);
edofVec = reshape(2 * nodenrs(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] -2 -1], nelx * nely, 1);
edofMat = [edofMat;
        267, 268, 275, 276, 273, 274, 265, 266;
        269, 270, 277, 278, 275, 276, 267, 268;
        271, 272, 279, 280, 277, 278, 269, 270;
        241, 242, 249, 250, 279, 280, 271, 272];
iK = reshape(kron(edofMat, ones(8, 1))', 64 * nel, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * nel, 1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = zeros(280, 1);
F([257 259 261 263]) = [-1 -2 -2 -1];
F([248 256]) = [-1/2 -1/2];
U = zeros(280, 1);
fixeddofs = [1:8, 266, 274];
alldofs = 1:280;
freedofs = setdiff(alldofs, fixeddofs);
%% PREPARE FILTER
iH = ones(nelx * nely * (2 * (ceil(rmin) - 1) + 1)^2, 1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;

%% INITIALIZE ITERATION
x = ones(nel, 1);
loop = 0;

dvol = 1;
stage = 1;
totalLoop = 0;
femAnalysis = 0;
epsilon = 1e-4;

while stage < 3
    loop = loop + 1;
    Lower = 0;
    Upper = 1e9;

    if stage == 1
        vol = floor(volfrac);
        vol = vol - dvol;
        vol = max(vol, volfrac0);
        volfrac = vol;
    end

    if stage == 2
        vol = floor(volfrac0);
        epsilon = 1e-4;
    end

    %         subplot(2, 1, 1);
    %     colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

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
%     x([2 5]) = 0;
    xPhys = x;
    xPhys(x < 1e-3) = Emin;
    sK = reshape(KE(:) * (xPhys(:)'), 64 * nel, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    ce = sum((U(edofMat) * KE) .* U(edofMat), 2);
    ce(x < 1e-3) = 3 * ce(x < 1e-3);
    c = F' * U/2;
    disp(c);
    %         ce(:) = ce .* x;
    %     ce(:) = H * (ce(:) ./ Hs);
    femAnalysis = femAnalysis + 1;

    tOpt1 = tic;
    [xResult, cost, exitFlag] = gbdMasterCut(reshape(x, [], 1), c, reshape(ce, 1, []), [], [], [], vol);
    optTime = optTime + toc(tOpt1);

    x = reshape(xResult, size(x, 1), size(x, 2));
    xOptimal = x;

    %         subplot(2, 1, 2);
    %     colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

    while (1)
        innerLoop = innerLoop + 1;
        disp((Upper - Lower) / Upper);

        % primal problem
        xPhys = x;
        xPhys(x < 1e-3) = Emin;
        sK = reshape(KE(:) * (xPhys(:)'), 64 * nel, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        ce = sum((U(edofMat) * KE) .* U(edofMat), 2);
        ce(x < 1e-3) = 3 * ce(x < 1e-3);
        c = F' * U/2;
        %                 ce(:) = ce .* x;
        %         ce(:) = H * (ce(:) ./ Hs);
        femAnalysis = femAnalysis + 1;

        if c < Upper
            xOptimal = x;
            Upper = c;
        end

        % check feasibility
        if norm(U) > 1e9
            U(freedofs) = feasibilityCut(K(freedofs, freedofs), F(freedofs), Upper);
            ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
            c = sum(sum((Emin + xPhys * (E0 - Emin)) .* ce));

            xFeasible = [xFeasible reshape(x, [], 1)];
            ceFeasible = [ceFeasible; reshape(ce, 1, [])];
            cFeasible = [cFeasible; c];
        else
            xTarget = [xTarget reshape(x, [], 1)];
            ceTarget = [ceTarget; reshape(ce, 1, [])];
            cTarget = [cTarget; c];
            index = [];

            for i = 1:length(cTarget)

                if (cTarget(i) <= c)
                    index = [index; i];
                end

            end

        end

        % master problem
        tOpt1 = tic;
        [xResult, cost, ~] = gbdMasterCut(xTarget(:, index), cTarget(index), ceTarget(index, :), xFeasible, cFeasible, ceFeasible, vol);
        tOp2 = toc(tOpt1);
        optTime = optTime + tOp2;

        x = xResult;

%         xPhys = x;
%         xPhys(x < 1e-3) = Emin;
%         sK = reshape(KE(:) * (Emin + xPhys(:)' * (E0 - Emin)), 64 * nel, 1);
%         K = sparse(iK, jK, sK); K = (K + K') / 2;
%         U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
%         cost = F' * U;

        %         subplot(2, 1, 1);
        %         colormap(gray); imagesc(1 - xOptimal); caxis([0 1]); axis equal; axis off; drawnow;
        %         subplot(2, 1, 2);
        xOutput = reshape(x(1:96), nely, nelx);
        colormap(gray); imagesc(1 - xOutput); caxis([0 1]); axis equal; axis off; drawnow;

        if cost > Upper || (Upper - cost) / Upper < epsilon
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f, Gap.:%5.3f%%\n', loop, Upper, sum(x(:)), (Upper - cost) / Upper * 100);
            break;
        end

        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f, Gap.:%5.3f%%\n', loop, Upper, sum(x(:)), (Upper - cost) / Upper * 100);
    end

    x = xOptimal;

    totalLoop = totalLoop + innerLoop;

    if stage == 1 && volfrac <= volfrac0
        stage = 2;
    elseif stage == 2
        stage = 3;
    end

end

disp(1 - x(97:end));

totalTime = toc(t1);
disp(['total time: ', num2str(totalTime), 's']);
disp(['optimization time: ', num2str(optTime), 's']);

% exit()

function [x, objFunc, exitFlag] = gbdMasterCut(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    n = size(y, 2);
    m = size(yFeasible, 2);

    if n > 1
        l = size(y, 1);
        f = zeros(1, l + 1);
        f(1) = 1;
        lb = zeros(1, l + 1);
        lb(1) = -inf;
        lb([95 97 98:end]) = 1;
        ub = ones(1, l + 1);
        ub(1) = inf;

        A = zeros(n + m, l + 1);
        b = zeros(n + m, 1);

        for i = 1:n
            A(i, 1) = -1;
            A(i, 2:end) = -weight(i, :);
            b(i) = -obj(i) - weight(i, :) * y(:, i);
        end

        for i = 1:m
            A(i + n, 2:end) = weightFeasible(i, :);
            b(i + n) = objFeasible(i);
        end

        intcon = 1:l;
        intcon = intcon + 1;

        Aeq = ones(1, l + 1);
        Aeq(1, 1) = 0;
        beq = vol;

        %     options = optimoptions('intlinprog','IntegerPreprocess','none');
        [x, objFunc, exitFlag, ~] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);

        if exitFlag ~= 1
            x = y;
        else
            x = x(2:end);
        end

    else
        l = size(y, 1);
        xOptimal = y(:, 1);
        objOptimal = inf;

        f = -weight(1, :);
        lb = zeros(1, l);
        ub = ones(1, l);
        lb([94 96 97:end]) = 1;

        intcon = 1:l;

        Aeq = ones(1, l);
        beq = vol;

        %     options = optimoptions('intlinprog','IntegerPreprocess','none');
        [x, ~, exitFlag, ~] = intlinprog(f, intcon, [], [], Aeq, beq, lb, ub);

        if exitFlag == 1
            objFunc = -inf;

            for j = 1:n
                objSelect = obj(j) - weight(j, :) * (x - y(:, j));

                if objSelect > objFunc
                    objFunc = objSelect;
                end

            end

            if objFunc < objOptimal
                xOptimal = x;
                objOptimal = objFunc;
            end

        end

        x = xOptimal;
        objFunc = objOptimal;
    end

end

function [x, objFunc, exitFlag] = gbdMasterCutPython(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    n = size(y, 2);
    m = size(yFeasible, 2);

    if n > 1
        l = size(y, 1);
        f = zeros(1, l + 1);
        f(1) = 1;
        lb = zeros(1, l + 1);
        lb(1) = -inf;
        ub = ones(1, l + 1);
        ub(1) = inf;

        A = zeros(n + m, l + 1);
        b = zeros(n + m, 1);

        for i = 1:n
            A(i, 1) = -1;
            A(i, 2:end) = -weight(i, :);
            b(i) = -obj(i) - weight(i, :) * y(:, i);
        end

        for i = 1:m
            A(i + n, 2:end) = weightFeasible(i, :);
            b(i + n) = objFeasible(i);
        end

        intcon = 1:l;
        intcon = intcon + 1;

        Aeq = ones(1, l + 1);
        Aeq(1, 1) = 0;
        beq = vol;

        save('cons_mip.mat', 'weight', 'b', 'vol');
        system('python3 cons_mip.py')
        load('cons_mip_result.mat');
        exitFlag = 1;

        objFunc = inf;

        %     options = optimoptions('intlinprog','IntegerPreprocess','none');
        [xCompare, objCompare, ~, ~] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);
        disp(norm(xCompare(2:end) - x, 1));

        for j = 1:n
            objSelect = obj(j) - weight(j, :) * (x - y(:, j));

            if objSelect < objFunc
                objFunc = objSelect;
            end

        end

        disp([objCompare, objFunc]);

        if exitFlag ~= 1
            x = y;
        else
            x = x;
        end

    else
        l = size(y, 1);
        xOptimal = y(:, 1);
        objOptimal = inf;

        f = -weight(1, :);
        lb = zeros(1, l);
        ub = ones(1, l);

        intcon = 1:l;

        Aeq = ones(1, l);
        beq = vol;

        save('mip.mat', 'f', 'vol', 'obj');
        system('python3 lin_mip.py')
        load('mip_result.mat');
        exitFlag = 1;

        %     options = optimoptions('intlinprog','IntegerPreprocess','none');
        % [x, ~, exitFlag, ~] = intlinprog(f, intcon, Aeq, beq, [], [], lb, ub);

        objFunc = -inf;

        objFunc = obj(1) - weight(1, :) * (x - y(:, 1));

        if objFunc < objOptimal
            xOptimal = x;
            objOptimal = objFunc;
        end

        x = xOptimal;
        objFunc = objOptimal;
    end

end

function lambda = feasibilityCut(K, f, UB)
    n = length(f);

    A = -f';
    b = -UB;

    lambda = fmincon(@(x)(0), zeros(n, 1), A, b, [], [], [], [], @con);

    function [c, ceq] = con(x)
        c = x' * K * x - f' * x;
        ceq = [];
    end

end
