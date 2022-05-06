%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function top_thermal_vol(nelx, nely, volfrac0, penal, rmin)
    close all;
    h1 = figure(1);
    set(h1, 'position', [10, 200, 600, 600])
    h2 = figure(2);
    set(h2, 'position', [700, 200, 600, 600])
    % INITIALIZE
    volfrac = 0.9;
    x(1:nely, 1:nelx) = volfrac;
    while volfrac >= volfrac0
    %     load('x_thermal.mat');
    loop = 0;
    change = 1.;

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
    F = sparse(1:(nelx + 1) * (nely + 1), 1, 0.01, (nelx + 1) * (nely + 1), 1);
    U = zeros((nely + 1) * (nelx + 1), 1);
    fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
    % fixeddofs = 1:nely+1;
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

    % START ITERATION
    while change > 0.01
        loop = loop + 1;
        xold = x;
        % FE-ANALYSIS
        xPhys = x;
        sK = reshape(KE(:) * (kmin + xPhys(:)'.^penal * (k0 - kmin)), 16 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
        c = sum(sum((kmin + xPhys.^penal * (k0 - kmin)) .* ce));
        dc = -penal * (k0 - kmin) * xPhys.^(penal - 1) .* ce;
%         [dc] = check(nelx, nely, rmin, x, dc);
        dc(:) = H * (dc(:) ./ Hs);
        % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
        [x] = OC(nelx, nely, x, volfrac, dc);
        % PRINT RESULTS
        change = max(max(abs(x - xold)));
        disp([' It.: ' sprintf('%4i', loop) ' Obj.: ' sprintf('%10.4f', c) ...
                ' Vol.: ' sprintf('%6.3f', sum(sum(x)) / (nelx * nely)) ...
                ' ch.: ' sprintf('%6.3f', change)])
        % PLOT DENSITIES
        figure(1);
        colormap(gray); imagesc(-x); axis equal; axis tight; axis off;
        figure(2);
%         xMesh = linspace(0, 1, nelx + 1);
%         yMesh = linspace(0, 1, nely + 1);
%         [xMesh, yMesh] = meshgrid(xMesh, yMesh);
        pcolor(reshape(U, nely + 1, nelx + 1));
        axis equal;
        view(0, 90);
        colorbar;

        pause(1e-6);
    end

%     sum(x(:));
% 
%         x(x>0.5) = 1;
%         x(x<0.5) = 0;
%         disp(sum(x(:)));
% 
%         [U] = FE(nelx, nely, x, penal);
%         figure;
%         xMesh = linspace(0, 1, nelx+1);
%         yMesh = linspace(0, 1, nely+1);
%         [xMesh, yMesh] = meshgrid(xMesh, yMesh);
%         mesh(xMesh, yMesh, reshape(U, nely+1, nelx+1));
%         view(0, 0);
%         colorbar;
%         % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
%         [KE] = lk;
%         c = 0.;
%     
%         for ely = 1:nely
%     
%             for elx = 1:nelx
%                 n1 = (nely + 1) * (elx - 1) + ely;
%                 n2 = (nely + 1) * elx +ely;
%                 Ue = U([n1; n2; n2 + 1; n1 + 1], 1);
%                 c = c + (0.0001 + 0.9999 * x(ely, elx)^penal) * Ue' * KE * Ue;
%                 dc(ely, elx) = -0.9999 * penal * x(ely, elx)^(penal - 1) * Ue' * KE * Ue;
%             end
%     
%         end
%         disp(c);
%         save('x_thermal', 'x');
    %     colormap(gray); imagesc(-x); axis equal; axis tight; axis off; pause(1e-6);
        volfrac = volfrac - 0.1;
        x = x / (volfrac + 0.1) * volfrac;
    end

    x(x>0.5) = 1;
    x(x<0.5) = 1e-3;
    xPhys = x;
    sK = reshape(KE(:) * (kmin + xPhys(:)'.^penal * (k0 - kmin)), 16 * nelx * nely, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
    c = sum(sum((kmin + xPhys.^penal * (k0 - kmin)) .* ce));
    disp(c);
    disp(sum(x(:)));
    figure(1);
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off;
end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew] = OC(nelx, nely, x, volfrac, dc)
    l1 = 0; l2 = 100000; move = 0.2;

    while (l2 - l1 > 1e-6)
        lmid = 0.5 * (l2 + l1);
        xnew = max(0.001, max(x - move, min(1., min(x + move, x .* sqrt(-dc ./ lmid)))));

        if sum(sum(xnew)) - volfrac * nelx * nely > 0
            l1 = lmid;
        else
            l2 = lmid;
        end

    end

end

%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U] = FE(nelx, nely, x, penal)
    [KE] = lk;
    K = sparse((nelx + 1) * (nely + 1), (nelx + 1) * (nely + 1));
    F = sparse((nely + 1) * (nelx + 1), 1); U = zeros((nely + 1) * (nelx + 1), 1);

    k = 1.0;
    kmin = 0.9;

    for elx = 1:nelx

        for ely = 1:nely
            n1 = (nely + 1) * (elx - 1) + ely;
            n2 = (nely + 1) * elx +ely;
            edof = [n1; n2; n2 + 1; n1 + 1];
            K(edof, edof) = K(edof, edof) + ((k - kmin) * x(ely, elx)^penal + kmin) * KE;
        end

    end

    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F(1:(nely + 1) * (nelx + 1), 1) = 1;
    %     F(ceil((nely+1)*(nelx+1)/2), 1) = -1.0;
    fixeddofs = nely / 2 + 1 - floor(nely / 20):nely / 2 + 1 + floor(nely / 20);
    alldofs = 1:(nely + 1) * (nelx + 1);
    freedofs = setdiff(alldofs, fixeddofs);
    % SOLVING
    U(freedofs, :) = K(freedofs, freedofs) \ F(freedofs, :);
    U(fixeddofs, :) = 0;
end

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE] = lk
    KE = [2/3 -1/6 -1/3 -1/6
        -1/6 2/3 -1/6 -1/3
        -1/3 -1/6 2/3 -1/6
        -1/6 -1/3 -1/6 2/3];
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
