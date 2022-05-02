clear all;
close all;
clc;

nelx = 40;
nely = 40;
nelz = 5;
rmin = 1.4;
volfrac0 = 0.3;
volfrac = 1.0;
%% MATERIAL PROPERTIES
k0 = 1; % Good thermal conductivity
kmin = 1e-3; % Poor thermal conductivity
% USER-DEFINED SUPPORT FIXED DOFs
il = nelx / 2 - nelx / 20:nelx / 2 + nelx / 20; jl = nely; kl = 0:nelz;
fixedxy = il * (nely + 1) + (nely + 1 - jl);
fixednid = repmat(fixedxy', size(kl)) + ...
    repmat(kl * (nelx + 1) * (nely + 1), size(fixedxy, 2), 1);
fixeddof = reshape(fixednid, [], 1);
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx * nely * nelz;
ndof = (nelx + 1) * (nely + 1) * (nelz + 1);
F = sparse(1:ndof, 1, -0.01, ndof, 1);
U = zeros(ndof, 1);
freedofs = setdiff(1:ndof, fixeddof);
KE = lk_H8(k0);
nodegrd = reshape(1:(nely + 1) * (nelx + 1), nely + 1, nelx + 1);
nodeids = reshape(nodegrd(1:end - 1, 1:end - 1), nely * nelx, 1);
nodeidz = 0:(nely + 1) * (nelx + 1):(nelz - 1) * (nely + 1) * (nelx + 1);
nodeids = repmat(nodeids, size(nodeidz)) + repmat(nodeidz, size(nodeids));
edofVec = nodeids(:) + 1;
edofMat = repmat(edofVec, 1, 8) + ...
    repmat([0 nely + [1 0] -1 ...
                (nely + 1) * (nelx + 1) + [0 nely + [1 0] -1]], nele, 1);
iK = reshape(kron(edofMat, ones(8, 1))', 8 * 8 * nele, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 8 * 8 * nele, 1);
% PREPARE FILTER
iH = ones(nele * (2 * (ceil(rmin) - 1) + 1)^2, 1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;

for k1 = 1:nelz

    for i1 = 1:nelx

        for j1 = 1:nely
            e1 = (k1 - 1) * nelx * nely + (i1 - 1) * nely + j1;

            for k2 = max(k1 - (ceil(rmin) - 1), 1):min(k1 + (ceil(rmin) - 1), nelz)

                for i2 = max(i1 - (ceil(rmin) - 1), 1):min(i1 + (ceil(rmin) - 1), nelx)

                    for j2 = max(j1 - (ceil(rmin) - 1), 1):min(j1 + (ceil(rmin) - 1), nely)
                        e2 = (k2 - 1) * nelx * nely + (i2 - 1) * nely + j2;
                        k = k + 1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0, rmin - sqrt((i1 - i2)^2 + (j1 - j2)^2 + (k1 - k2)^2));
                    end

                end

            end

        end

    end

end

H = sparse(iH, jH, sH);
Hs = sum(H, 2);
%% INITIALIZE ITERATION
x = ones(nely, nelx, nelz);
loop = 0;

dvol = 200;
stage = 1;
totalLoop = 0;
epsilon = 1e-2;

numFeasibleCut = 0;

compliance = [];
volume = [];

while stage < 3
    loop = loop + 1;
    Lower = 0;
    Upper = 1e9;

    if stage == 1
        vol = floor(volfrac * nelx * nely * nelz);
        vol = vol - dvol;
        volfrac = vol / (nelx * nely * nelz);
    end

    if stage == 2
        vol = floor(volfrac0 * nelx * nely * nelz);
        epsilon = 1e-5;
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
    xPhys = x;
    sK = reshape(KE(:) * (kmin + (1 - kmin) * xPhys(:)'), 8 * 8 * nele, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    U(freedofs, :) = K(freedofs, freedofs) \ F(freedofs, :);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), [nely, nelx, nelz]);
    c = sum(sum(sum((kmin + (1 - kmin) * xPhys) .* ce)));
    ce(:) = ce .* x;
    ce(:) = H * (ce(:) ./ Hs);

    [xResult, cost, exitFlag] = gbdMasterCut(reshape(x, [], 1), c, reshape(ce, 1, []), [], [], [], vol);
    %     [xResult, cost, exitFlag] = gbdMasterCutRelaxed(reshape(x, [], 1), c, reshape(ce, 1, []), [], [], [], vol);
    % [xResult, cost, exitFlag] = gbdMasterCutQuantum(reshape(x, [], 1), c, reshape(ce, 1, []), [], [], [], vol);

    x = reshape(xResult, size(x, 1), size(x, 2), size(x, 3));
    xOptimal = x;

    while (1)
        innerLoop = innerLoop + 1;

        % primal problem
        xPhys = x;
        sK = reshape(KE(:) * (kmin + (1 - kmin) * xPhys(:)'), 8 * 8 * nele, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs, :) = K(freedofs, freedofs) \ F(freedofs, :);
        % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), [nely, nelx, nelz]);
        c = sum(sum(sum((kmin + (1 - kmin) * xPhys) .* ce)));
        ce(:) = ce .* x;
        ce(:) = H * (ce(:) ./ Hs);

        if c < Upper
            xOptimal = x;
            Upper = c;
        end

        compliance = [compliance; Upper];
        volume = [volume; vol];

        % check feasibility
        if norm(U) > 1e9
            numFeasibleCut = numFeasibleCut + 1;
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
        % if stage == 1
        [xResult, cost, exitFlag] = gbdMasterCut(xTarget(:, index), cTarget(index), ceTarget(index, :), xFeasible, cFeasible, ceFeasible, vol);
        %         [xResult, cost, exitFlag] = gbdMasterCutRelaxed(xTarget(:, index), cTarget(index), ceTarget(index, :), xFeasible, cFeasible, ceFeasible, vol);
        % else
        %     [xResult, cost, exitFlag] = gbdMasterCutQuantum(xTarget(:, index), cTarget(index), ceTarget(index, :), xFeasible, cFeasible, ceFeasible, vol);
        % end

        x = reshape(xResult, size(x, 1), size(x, 2), size(x, 3));

        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f, Gap.:%5.3f%%\n', loop, Upper, sum(x(:)), (Upper - cost) / Upper * 100);

        if cost > Upper || (Upper - cost) / Upper < epsilon
            break;
        end

    end

    x = xOptimal;

    totalLoop = totalLoop + innerLoop;

    if stage == 1 && volfrac <= volfrac0
        stage = 2;
    elseif stage == 2
        stage = 3;
    end

end

display_3D(xPhys);

disp(numFeasibleCut)
disp(Upper)
disp(totalLoop)

figure;
hold on;
yyaxis left
plot(compliance, 'b.-');
yyaxis right
plot(volume ./ (nelx * nely * nelz), 'r.-');

function [x, objFunc, exitFlag] = gbdMasterCut(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    n = size(y, 2);
    m = size(yFeasible, 2);
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

    %     options = optimoptions('intlinprog','IntegerPreprocess','none');
    [x, objFunc, exitFlag, ~] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);

    if exitFlag ~= 1
        x = y;
    else
        x = x(2:end);
    end

end

function [KE] = lk_H8(k)
    A1 = 4 * eye(2); A2 = -eye(2);
    A3 = fliplr(A2); A4 = -ones(2);
    KE1 = [A1 A2; A2 A1];
    KE2 = [A3 A4; A4 A3];
    KE = 1/12 * k * [KE1 KE2; KE2 KE1];
end

function [x, objFunc, exitFlag] = gbdMasterCutRelaxed(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    n = size(y, 2);
    m = size(yFeasible, 2);
    l = size(y, 1);

    xOptimal = y(:, 1);
    objOptimal = inf;

    for i = 1:n
        f = -weight(i, :);
        lb = zeros(1, l);
        ub = ones(1, l);

        intcon = 1:l;

        Aeq = ones(1, l);
        beq = vol;

        %     options = optimoptions('intlinprog','IntegerPreprocess','none');
        [x, ~, exitFlag, ~] = intlinprog(f, intcon, Aeq, beq, [], [], lb, ub);

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

    end

    x = xOptimal;
    objFunc = objOptimal;
end

function [x, objFunc, exitFlag] = gbdMasterCutQuantum(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    exitFlag = 1;

    % reduce the amount of variables
    weightSparse = sparse(weight);
    [~, nonzero] = find(weightSparse);
    nonzero = unique(sort(nonzero));

    weightReduced = weight(:, nonzero);
    objReduced = obj;
    yReduced = y(nonzero, :);

    save('quantum.mat', 'weightReduced', 'objReduced', 'yReduced', 'vol');
    system('python3 cut_solver.py');
    load('result.mat');

    x = zeros(size(y, 1), 1);
    x(nonzero) = res;

    [xCompare, objFuncCompare, ~] = gbdMasterCut(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol);
    fprintf("norm difference: %f\n", norm(x - xCompare, 1));
    disp([objFunc, objFuncCompare]);
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

function display_3D(rho)
    [nely, nelx, nelz] = size(rho);
    hx = 1; hy = 1; hz = 1; % User-defined unit element size
    face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
    set(gcf, 'Name', 'ISO display', 'NumberTitle', 'off');

    for k = 1:nelz
        z = (k - 1) * hz;

        for i = 1:nelx
            x = (i - 1) * hx;

            for j = 1:nely
                y = nely * hy - (j - 1) * hy;

                if (rho(j, i, k) > 0.5) % User-defined display density threshold
                    vert = [x y z; x y - hx z; x + hx y - hx z; x + hx y z; x y z + hx; x y - hx z + hx; x + hx y - hx z + hx; x + hx y z + hx];
                    vert(:, [2 3]) = vert(:, [3 2]); vert(:, 2, :) = -vert(:, 2, :);
                    patch('Faces', face, 'Vertices', vert, 'FaceColor', [0.2 + 0.8 * (1 - rho(j, i, k)), 0.2 + 0.8 * (1 - rho(j, i, k)), 0.2 + 0.8 * (1 - rho(j, i, k))]);
                    hold on;
                end

            end

        end

    end

    axis equal; axis tight; axis off; box on; view([30, 30]); pause(1e-6);
end
