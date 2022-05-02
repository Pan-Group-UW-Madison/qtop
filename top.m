%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function top(nelx, nely, volfrac, penal, rmin)
    h1 = figure(1);
    set(h1, 'position', [10, 200, 600, 600])
    h2 = figure(2);
    set(h2, 'position', [700, 200, 600, 600])
    % INITIALIZE
    x(1:nely, 1:nelx) = volfrac;
%     load('x_thermal.mat');
    loop = 0;
    change = 1.;
    % START ITERATION
    while change > 0.01
        loop = loop + 1;
        xold = x;
        % FE-ANALYSIS
        [U] = FE(nelx, nely, x, penal);
        % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        [KE] = lk;
        c = 0.;

        for ely = 1:nely

            for elx = 1:nelx
                n1 = (nely + 1) * (elx - 1) + ely;
                n2 = (nely + 1) * elx +ely;
                Ue = U([n1; n2; n2 + 1; n1 + 1], 1);
                c = c + (x(ely, elx)^penal) * Ue' * KE * Ue;
                dc(ely, elx) = -0.9999 * penal * x(ely, elx)^(penal - 1) * Ue' * KE * Ue;
            end

        end

        % FILTERING OF SENSITIVITIES
        [dc] = check(nelx, nely, rmin, x, dc);
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
%         xMesh = linspace(0, 1, nelx+1);
%         yMesh = linspace(0, 1, nely+1);
%         [xMesh, yMesh] = meshgrid(xMesh, yMesh);
        pcolor(reshape(U, nely+1, nelx+1));
        axis equal;
%         view(0, 90);
        colorbar;
%         caxis([0 2])
        
        pause(1e-6);
    end
    
    sum(x(:));

    
%     x(x>0.5) = 1;
%     x(x<0.5) = 0;
%     disp(sum(x(:)));
    
%     [U] = FE(nelx, nely, x, penal);
%     figure;
%     xMesh = linspace(0, 1, nelx+1);
%     yMesh = linspace(0, 1, nely+1);
%     [xMesh, yMesh] = meshgrid(xMesh, yMesh);
%     mesh(xMesh, yMesh, reshape(U, nely+1, nelx+1));
%     view(0, 0);
%     colorbar;
%     % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
%     [KE] = lk;
%     c = 0.;
% 
%     for ely = 1:nely
% 
%         for elx = 1:nelx
%             n1 = (nely + 1) * (elx - 1) + ely;
%             n2 = (nely + 1) * elx +ely;
%             Ue = U([n1; n2; n2 + 1; n1 + 1], 1);
%             c = c + (0.0001 + 0.9999 * x(ely, elx)^penal) * Ue' * KE * Ue;
%             dc(ely, elx) = -0.9999 * penal * x(ely, elx)^(penal - 1) * Ue' * KE * Ue;
%         end
% 
%     end
%     disp(c);
%     save('x_thermal', 'x');
%     colormap(gray); imagesc(-x); axis equal; axis tight; axis off; pause(1e-6);
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

    for elx = 1:nelx

        for ely = 1:nely
            n1 = (nely + 1) * (elx - 1) + ely;
            n2 = (nely + 1) * elx +ely;
            edof = [n1; n2; n2 + 1; n1 + 1];
            K(edof, edof) = K(edof, edof) + (x(ely, elx)^penal + 1e-3) * KE;
        end

    end

    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F(1:(nely + 1) * (nelx + 1), 1) = 0.01;
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
