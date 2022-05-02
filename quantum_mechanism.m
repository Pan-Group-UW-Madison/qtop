clear all;
close all;
clc;

figure;
set(gcf, 'position', [200, 200, 400, 200])

nelx = 40;
nely = 20;
rmin = 2;
volfrac0 = 0.3;
volfrac = 1.0;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 0.5;
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

din = 1;
dout = 2 * nelx * (nely + 1) + 1;
F = sparse(2 * (nely + 1) * (nelx + 1), 2);
U = zeros(2 * (nely + 1) * (nelx + 1), 2);
F(din, 1) = 1;
F(dout, 2) = -1;
fixeddofs = union([2:2 * (nely + 1):2 * (nelx + 1) * (nely + 1)], [2 * (nely + 1):-1:2 * (nely + 1) - 3]);

alldofs = [1:2 * (nely + 1) * (nelx + 1)];
freedofs = setdiff(alldofs, fixeddofs);
%% PREPARE FILTER
iH = ones(nelx * nely * (2 * (ceil(rmin) - 1) + 1)^2, 1);
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
x = ones(nely, nelx);
loop = 0;

dvol = 50;
stage = 1;
totalLoop = 0;
epsilon = 1e-3;

% stage = 2;
% load('xm');
% vol = floor(volfrac0*nelx*nely);

while stage < 3
    loop = loop + 1;
    Upper = 1e9;

    if stage == 1
        vol = floor(volfrac * nelx * nely);
        vol = vol - dvol;
        volfrac = vol / (nelx * nely);
    end

    if stage == 2
        epsilon = 1e-3;
    end

    %     vol = 1350;

    %     load('x.mat');
    %     x(xPhys>0.5) = 1;
    %     x(xPhys<0.5) = 0;

    %     subplot(2, 1, 1);
    colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

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
    %     xPhys = x;
    %     sK = reshape(KE(:)*(Emin+xPhys(:)'*(E0-Emin)),64*nelx*nely,1);
    %     K = sparse(iK,jK,sK); K = (K+K')/2;
    %     K(din, din) = K(din,din) +0.1;
    %     K(dout, dout) = K(dout,dout) +0.1;
    %     U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %     U(fixeddofs,:)= 0;
    %     ce=zeros(nely, nelx);
    %     for ely = 1:nely
    %         for elx = 1:nelx
    %           n1 = (nely+1)*(elx-1)+ely;
    %           n2 = (nely+1)* elx   +ely;
    %           Ue1 = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
    %           Ue2 = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],2);
    %           ce(ely,elx) = x(ely,elx)*Ue1'*KE*Ue2;
    %         end
    %     end
    %     c = U(dout, 2);
    %     ce(:) = ce.*x;
    %     ce(:) = H*(ce(:)./Hs);

    xPhys = x;
    xPhys(xPhys < 1e-3) = Emin;
    [U, c] = FE(nelx, nely, xPhys, 1.0);
    [KE] = lk;
    ce = zeros(nely, nelx);

    for ely = 1:nely

        for elx = 1:nelx
            n1 = (nely + 1) * (elx - 1) + ely;
            n2 = (nely + 1) * elx +ely;
            Ue1 = U([2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2; 2 * n2 + 1; 2 * n2 + 2; 2 * n1 + 1; 2 * n1 + 2], 1);
            Ue2 = U([2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2; 2 * n2 + 1; 2 * n2 + 2; 2 * n1 + 1; 2 * n1 + 2], 2);
            ce(ely, elx) =- Ue1' * KE * Ue2;
        end

    end

    %     ce(:) = 3 * ce .* x.^2;
    %     ce(x<1e-3) = ce(x<1e-3)*Emin^2;
    ce(:) = H * (ce(:) ./ Hs);

    [xResult, cost, exitFlag] = gbdMasterCut(reshape(x, [], 1), c, reshape(ce, 1, []), [], [], [], vol);

    x = reshape(xResult, size(x, 1), size(x, 2));
    xOptimal = x;

    %     subplot(2, 1, 2);
    colormap(gray); imagesc(1 - x); caxis([0 1]); axis equal; axis off; drawnow;

    while (1)
        innerLoop = innerLoop + 1;

        % primal problem
        %         xPhys = x;
        %         sK = reshape(KE(:)*(Emin+xPhys(:)'*(E0-Emin)),64*nelx*nely,1);
        %         K = sparse(iK,jK,sK); K = (K+K')/2;
        %         K(din, din) = K(din,din) +0.1;
        %         K(dout, dout) = K(dout,dout) +0.1;
        %         U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
        %         U(fixeddofs,:)= 0;
        %         ce=zeros(nely, nelx);
        %         for ely = 1:nely
        %             for elx = 1:nelx
        %               n1 = (nely+1)*(elx-1)+ely;
        %               n2 = (nely+1)* elx   +ely;
        %               Ue1 = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
        %               Ue2 = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],2);
        %               ce(ely,elx) = x(ely,elx)*Ue1'*KE*Ue2;
        %             end
        %         end
        %         c = U(dout, 2);
        %         ce(:) = ce.*x;
        %         ce(:) = H*(ce(:)./Hs);

        xPhys = x;
        xPhys(xPhys < 1e-3) = Emin;
        [U, c] = FE(nelx, nely, xPhys, 1.0);
        [KE] = lk;
        ce = zeros(nely, nelx);

        for ely = 1:nely

            for elx = 1:nelx
                n1 = (nely + 1) * (elx - 1) + ely;
                n2 = (nely + 1) * elx +ely;
                Ue1 = U([2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2; 2 * n2 + 1; 2 * n2 + 2; 2 * n1 + 1; 2 * n1 + 2], 1);
                Ue2 = U([2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2; 2 * n2 + 1; 2 * n2 + 2; 2 * n1 + 1; 2 * n1 + 2], 2);
                ce(ely, elx) =- Ue1' * KE * Ue2;
            end

        end

        %         ce(x<1e-3) = ce(x<1e-3)*Emin^2;
        ce(:) = H * (ce(:) ./ Hs);

        if c < Upper
            xOptimal = x;
            Upper = c;
        end

        % check feasibility
        if norm(U(:, 1)) > 1e9
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
        [xResult, cost, exitFlag] = gbdMasterCut(xTarget(:, index), cTarget(index), ceTarget(index, :), xFeasible, cFeasible, ceFeasible, vol);

        if exitFlag == 1
            x = reshape(xResult, size(x, 1), size(x, 2));
        else
            break;
        end
        
%         xPhys = x;
%         xPhys(xPhys < 1e-3) = Emin;
%         [U, cost] = FE(nelx, nely, xPhys, 1.0);

        %         subplot(2, 1, 1);
        colormap(gray); imagesc(1 - xOptimal); caxis([0 1]); axis equal; axis off; drawnow;
        %         subplot(2, 1, 2);
        %         colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f, Gap.:%5.3f%%\n', loop, Upper, sum(x(:)), (Upper - cost) / abs(Upper) * 100);

        if cost > Upper || (Upper - cost) / abs(Upper) < epsilon
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

function [x, objFunc, exitFlag] = gbdMasterCut(y, obj, weight, yFeasible, objFeasible, weightFeasible, vol)
    n = size(y, 2);
    m = size(yFeasible, 2);
    l = size(y, 1);
    f = zeros(1, l + 1);
    f(1) = 1;
    lb = zeros(1, l + 1);
    lb(1) = -inf;
    lb(2) = 1;
    lb(1 + 20 * 39 + 1) = 1;
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

function [U, Uout] = FE(nelx, nely, x, penal)
    [KE] = lk;
    K = sparse(2 * (nelx + 1) * (nely + 1), 2 * (nelx + 1) * (nely + 1));
    F = sparse(2 * (nely + 1) * (nelx + 1), 2); U = sparse(2 * (nely + 1) * (nelx + 1), 2);

    for elx = 1:nelx

        for ely = 1:nely
            n1 = (nely + 1) * (elx - 1) + ely;
            n2 = (nely + 1) * elx +ely;
            edof = [2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2; 2 * n2 + 1; 2 * n2 + 2; 2 * n1 + 1; 2 * n1 + 2];
            K(edof, edof) = K(edof, edof) + (x(ely, elx)^penal + 1e-3) * KE;
        end

    end

    % DEFINE LOADS AND SUPPORTS (HALF FORCE INVERTER)
    din = 1;
    dout = 2 * nelx * (nely + 1) + 1;
    F(din, 1) = 1;
    F(dout, 2) = -1;
    K(din, din) = K(din, din) +0.1;
    K(dout, dout) = K(dout, dout) +0.1;
    fixeddofs = union([2:2 * (nely + 1):2 * (nelx + 1) * (nely + 1)], [2 * (nely + 1):-1:2 * (nely + 1) - 3]);

    alldofs = [1:2 * (nely + 1) * (nelx + 1)];
    freedofs = setdiff(alldofs, fixeddofs);
    % SOLVING
    U(freedofs, :) = K(freedofs, freedofs) \ F(freedofs, :);
    U(fixeddofs, :) = 0;
    Uout = full(U(dout, 1));
end

function [KE] = lk
    E = 1.;
    nu = 0.3;
    k = [1/2 - nu / 6 1/8 + nu / 8 -1/4 - nu / 12 -1/8 + 3 * nu / 8 ...
            -1/4 + nu / 12 -1/8 - nu / 8 nu / 6 1/8 - 3 * nu / 8];
    KE = E / (1 - nu^2) * [k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                        k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                        k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                        k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                        k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                        k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                        k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                        k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
