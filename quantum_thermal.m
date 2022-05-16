clear all;
close all;
clc;

figure;
set(gcf, 'position', [200, 200, 600, 600]);
optTime = 0;
t1 = tic;

numFEM = 0;
scale = 1e9;
nelx = 200;
nely = 200;
rmin = 1.5;
volfrac0 = 0.5;
volfrac = 1.0;
%% MATERIAL PROPERTIES
k0 = 1.0; % good thermal conductivity
kmin = 1e-3; % poor thermal conductivity
%% PREPARE FINITE ELEMENT ANALYSIS
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
F = sparse(1:(nelx + 1) * (nely + 1), 1, 1 / nelx / nely, (nelx + 1) * (nely + 1), 1);
U = zeros((nely + 1) * (nelx + 1), 1);
fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
x = ones(nely, nelx);
sK = reshape(KE(:) * (x(:)'), 16 * nelx * nely, 1);
K = sparse(iK, jK, sK); K = (K + K') / 2;
alldofs = 1:(nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);
%% PREPARE FILTER
iH = ones(nelx * nely * ((ceil(rmin) - 1) + 1)^2, 1);
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
                sH(k) = max(0, rmin - sqrt((i1 - i2)^2 + (j1 - j2)^2));
            end

        end

    end

end

H = sparse(iH, jH, sH);
Hs = sum(H, 2);
%% INITIALIZE ITERATION
loop = 0;
x = ones(nely, nelx);

dvol = 0.98;
stage = 1;
totalLoop = 0;
epsilon = 1e-7;

numFeasibleCut = 0;

compliance = [];
volume = [];

% stage = 2;
% load('x_thermal');
% 
% xPhys = x;
% 
% level = 0.5;
% levelHigh = 1;
% levelLow = kmin;
% while (true)
%     x = xPhys;
%     x(xPhys>=level) = 1;
%     x(xPhys<level) = kmin;
%     if (abs(sum(sum(x)) / (nely * nelx) - volfrac0) < (1/nely/nelx))
%         break;
%     else
%         if sum(sum(x)) / (nely * nelx) < volfrac0
%             levelHigh = level;
%             level = 0.5 * (levelLow + levelHigh);
%         else
%             levelLow = level;
%             level = 0.5 * (levelLow + levelHigh);
%         end
%     end
% end

% 
% x(x>=0.5) = 1;
% x(x<0.5) = 0;

% colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

while stage < 3
    loop = loop + 1;
    Lower = 0;
    Upper = 1e11;

    if stage == 1
        vol = floor(volfrac * nelx * nely);
        vol = floor(dvol * vol);
        if mod(vol, 2) == 1
            vol = vol - 1;
        end
%         vol = vol - dvol;
        vol = max(vol, volfrac0 * nelx * nely);
        volfrac = vol / (nelx * nely);
    end

    if stage == 2
        vol = floor(volfrac0 * nelx * nely);
%         kmin = 1e-3;
    end

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
    if stage == 1
        xPhys = x;
        xPhys(x < 1e-3) = kmin;
        xPhys(x > 0.5) = k0;

%         xPhys = x;
%         xPhys(:) = H * (xPhys(:) ./ Hs);
%         xPhys(xPhys > 1) = 1;
%         xPhys(xPhys < 1e-3) = kmin;

        sK = reshape(KE(:) * (xPhys(:)'), 16 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
%         
        ce(x < 0.5) = ce(x < 0.5) .* xPhys(x < 0.5).^2;
%         ce(:) = ce .* x;
%         ce(:) = H * (ce(:) ./ Hs);
        c = F' * U;
%         if stage == 2
%             
%         end
        ce = ce * scale;
        c = c * scale;
        numFEM = numFEM + 1;
    
        optTimer1 = tic;
        [xResult, cost, exitFlag] = gbdMasterCut(reshape(x, [], 1), c, reshape(ce, 1, []), vol);
        optTime = optTime + toc(optTimer1);
    
        x = reshape(xResult, size(x, 1), size(x, 2));
        xOptimal = x;
    end

%     colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

    while (1)
        innerLoop = innerLoop + 1;

        % primal problem
        xPhys = x;
        xPhys(x < 1e-3) = kmin;
        xPhys(x > 0.5) = k0;

%         xPhys = x;
%         xPhys(:) = H * (xPhys(:) ./ Hs);
%         xPhys(xPhys > 1) = 1;
%         xPhys(xPhys < 1e-3) = kmin;

        sK = reshape(KE(:) * (xPhys(:)'), 16 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
%         ce(:) = H * (ce(:) ./ Hs);
        ce(x < 0.5) = ce(x < 0.5) .* xPhys(x < 0.5).^2;
%         ce(:) = ce .* x;
%         ce(:) = H * (ce(:) ./ Hs);
        c = F' * U;
%         if stage == 2
%             ce(:) = H * (ce(:) ./ Hs);
%         end
        ce = ce * scale;
        numFEM = numFEM + 1;

        if c < Upper
            xOptimal = x;
            Upper = c;
        end
        c = c * scale;

        % check feasibility
        if norm(U) > 1e9
            numFeasibleCut = numFeasibleCut + 1;
            U(freedofs) = feasibilityCut(K(freedofs, freedofs), F(freedofs), Upper);
            ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
            c = sum(sum((Emin + xPhys * (k0 - kmin)) .* ce));

            xFeasible = [xFeasible reshape(x, [], 1)];
            ceFeasible = [ceFeasible; reshape(ce, 1, [])];
            cFeasible = [cFeasible; c];
        else
            xTarget = [xTarget reshape(x, [], 1)];
            ceTarget = [ceTarget; reshape(ce, 1, [])];
            cTarget = [cTarget; c];
            index = [];
%             index = 1:length(cTarget);

            for i = 1:length(cTarget)
                if (cTarget(i) <= c)
                    index = [index; i];
                end
            end

        end

        % master problem
        optTimer1 = tic;
        [xResult, cost, exitFlag] = gbdMasterCut(xTarget(:, index), cTarget(index), ceTarget(index, :), vol);
        cost = cost / scale;
        optTime = optTime + toc(optTimer1);

        x = reshape(xResult, size(x, 1), size(x, 2));

        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f, Gap.:%5.3f%%, index: %d\n', loop, Upper, sum(x(:)/nelx/nely), (Upper - cost) / Upper * 100, length(index));

%         figure(1);
%         colormap(gray); imagesc(1 - xOptimal); caxis([0 1]); axis equal; axis off; drawnow;
%         figure(2);
%         Ue = reshape(U, nely + 1, nelx + 1);
%         Ue = Ue(end:-1:1, :);
%         h = pcolor(Ue);
%         set(h, 'EdgeColor', 'none');
%         axis equal;
%         colorbar;
%         caxis([0 5]);
%         figure(3);
%         hold on;
%         yyaxis left
%         plot(compliance, 'b.-');
%         yyaxis right
%         plot(volume ./ (nelx * nely), 'r.-');

        if cost > Upper || (Upper - cost) / Upper < epsilon
            compliance = [compliance; Upper];
            volume = [volume; vol];
            break;
        end

        compliance = [compliance; Upper];
        volume = [volume; vol];

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

figure(1);
colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

figure(3);
hold on;
yyaxis left
plot(compliance, 'b.-');
yyaxis right
plot(volume ./ (nelx * nely), 'r.-');

figure(2);
Ue = reshape(U, nely + 1, nelx + 1);
Ue = Ue(end:-1:1, :);
h = pcolor(Ue);
set(h, 'EdgeColor', 'none');
axis equal;
colorbar;
caxis([0 2]);

% save('x_thermal', 'x');

function [x, objFunc, exitFlag] = gbdMasterCut(y, obj, weight, vol)
    n = size(y, 2);
    if n > 1
        l = size(y, 1);
        f = zeros(1, l + 1);
        f(1) = 1;
        lb = zeros(1, l + 1);
        lb(1) = 0;
        ub = ones(1, l + 1);
        ub(1) = inf;

        A = zeros(n + 1, l + 1);
        b = zeros(n + 1, 1);

        for i = 1:n
            A(i, 1) = -1;
            A(i, 2:end) = -weight(i, :);
            b(i) = -obj(i) - weight(i, :) * y(:, i);
        end

        intcon = 1:l;
        intcon = intcon + 1;

        A(end, :) = ones(1, l + 1);
        A(end, 1) = 0;
        b(end) = vol;

        [x, objFunc, exitFlag, ~] = intlinprog(f, intcon, A, b, [], [], lb, ub);

        x = x(2:end);

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

        x = xOptimal;
        objFunc = objOptimal;
    end

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

        A = ones(1, l);
        b = vol;

        %     options = optimoptions('intlinprog','IntegerPreprocess','none');
        [x, ~, exitFlag, ~] = intlinprog(f, intcon, [], [], A, b, lb, ub);

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

function [dcn] = check(nelx, nely, rmin, x, dc)
    dcn = zeros(nely, nelx);

    for i = 1:nelx

        for j = 1:nely
            sum = 0.0;

            for k = max(i - floor(rmin), 1):min(i + floor(rmin), nelx)

                for l = max(j - floor(rmin), 1):min(j + floor(rmin), nely)
                    fac = rmin - sqrt((i - k)^2 + (j - l)^2);
                    sum = sum + max(0, fac);
                    dcn(j, i) = dcn(j, i) + max(0, fac) * x(l, k) * dc(l, k);
                end

            end

            dcn(j, i) = dcn(j, i) / (x(j, i) * sum);
        end

    end

end
