function top_thermal_h(nelx, nely, volfrac, penal, rmin, ft)
    optTime = 0;
    t1 = tic;
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
    % fixeddofs = 1:nely+1;
    alldofs = 1:(nely + 1) * (nelx + 1);
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
    x = repmat(volfrac, nely, nelx);
    beta = 1;

    if ft == 1 || ft == 2
        xPhys = x;
    elseif ft == 3
        xTilde = x;
        xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta);
    end

    loopbeta = 0;
    loop = 0;
    change = 1;
    %% START ITERATION
    while change > 0.01 && loop < 100
        loopbeta = loopbeta + 1;
        loop = loop + 1;
        %% FE-ANALYSIS
        sK = reshape(KE(:) * (kmin + xPhys(:)'.^penal * (k0 - kmin)), 16 * nelx * nely, 1);
        K = sparse(iK, jK, sK); K = (K + K') / 2;
        U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
        c = sum(sum((kmin + xPhys.^penal * (k0 - kmin)) .* ce));
        dc = -penal * (k0 - kmin) * xPhys.^(penal - 1) .* ce;
        dv = ones(nely, nelx);
        %% FILTERING/MODIFICATION OF SENSITIVITIES
        if ft == 1
            dc(:) = H * (x(:) .* dc(:)) ./ Hs ./ max(1e-3, x(:));
        elseif ft == 2
            dc(:) = H * (dc(:) ./ Hs);
            dv(:) = H * (dv(:) ./ Hs);
        elseif ft == 3
            dx = beta * exp(-beta * xTilde) + exp(-beta);
            dc(:) = H * (dc(:) .* dx(:) ./ Hs);
            dv(:) = H * (dv(:) .* dx(:) ./ Hs);
        end

        %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e9; move = 0.2;
        optTimer1 = tic;

        while (l2 - l1) / (l1 + l2) > 1e-3
            lmid = 0.5 * (l2 + l1);
            xnew = max(0, max(x - move, min(1, min(x + move, x .* sqrt(-dc ./ dv / lmid)))));

            if ft == 1
                xPhys = xnew;
            elseif ft == 2
                xPhys(:) = (H * xnew(:)) ./ Hs;
            elseif ft == 3
                xTilde(:) = (H * xnew(:)) ./ Hs;
                xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta);
            end

            if sum(xPhys(:)) > volfrac * nelx * nely, l1 = lmid; else l2 = lmid; end
        end

        optTime = optTime + toc(optTimer1);
        change = max(abs(xnew(:) - x(:)));
        x = xnew;
        %% PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n', loop, c, ...
        mean(xPhys(:)), change);
        %% PLOT DENSITIES
%         colormap(gray); imagesc(1 - xPhys); caxis([0 1]); axis equal; axis off; drawnow;
        %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
        if ft == 3 && beta < 1024 && (loopbeta >= 50 || change <= 0.01)
            beta = 2 * beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n', beta);
        end

    end

    totalTime = toc(t1);
    disp(['total time: ', num2str(totalTime), 's']);
    disp(['optimization time: ', num2str(optTime), 's']);

    close all;
    grey_ratio = sum(sum((xPhys>=0.999)+(xPhys<=0.001))) / (nelx * nely);
    disp(['grey ratio: ', num2str(1 - grey_ratio)]);
    colormap(gray); imagesc(1 - xPhys); caxis([0 1]); axis equal; axis off; drawnow;

    x = xPhys;
    save('x_thermal', 'x');

    level = 0.5;
    levelHigh = 1;
    levelLow = kmin;
    while (true)
        x = xPhys;
        x(xPhys>=level) = 1;
        x(xPhys<level) = kmin;
        if (abs(sum(sum(x)) / (nely * nelx) - volfrac) < (1/nely/nelx))
            break;
        else
            if sum(sum(x)) / (nely * nelx) < volfrac
                levelHigh = level;
                level = 0.5 * (levelLow + levelHigh);
            else
                levelLow = level;
                level = 0.5 * (levelLow + levelHigh);
            end
        end
    end

    xPhys = x;

    sK = reshape(KE(:) * (kmin + xPhys(:)'.^penal * (k0 - kmin)), 16 * nelx * nely, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), nely, nelx);
    c = sum(sum((kmin + xPhys.^penal * (k0 - kmin)) .* ce));

    disp(c);
    disp(sum(sum(xPhys)) / (nelx * nely));
    figure;
    colormap(gray); imagesc(1 - xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    figure;
    h = pcolor(reshape(U, nelx + 1, nely + 1));
    set(h, 'EdgeColor', 'none');
    colorbar;
    axis equal;

%     save('x_thermal', 'x');
end