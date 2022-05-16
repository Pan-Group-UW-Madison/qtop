%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
%%%% Compliant Mechanism Design %%%%%
% topm(40,20,0.3,3.0, 1.2)
function topm(nelx,nely,volfrac,penal,rmin)
% INITIALIZE
x(1:nely,1:nelx) = volfrac; 

optTime = 0;
t1 = tic;

E = 1.; 
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

E0 = 1;
Emin = 1e-9;
F = zeros(2*(nely+1)*(nelx+1),2);
U = zeros(2*(nely+1)*(nelx+1),2);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([-2 -1 2*nely+[0 1 2 3] 0 1],nelx*nely,1);
iK = reshape(kron(edofMat, ones(8, 1))', 64 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * nelx * nely, 1);
% DEFINE LOADS AND SUPPORTS (HALF FORCE INVERTER)
din = 1;
dout = 2*nelx*(nely+1)+1;
F(din, 1) = 1;
F(dout, 2) = -1;
fixeddofs   = union([2:2*(nely+1):2*(nelx+1)*(nely+1)],[2*(nely+1):-1:2*(nely+1)-3]);
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);

loop = 0; 
change = 1.;
% START ITERATION
while change > 0.01 && loop < 500 
  loop = loop + 1;
  xold = x;
  
% FE-ANALYSIS
%     [U, c]=FE(nelx,nely,x,penal);  % displacement & cost function (output displcement)

    sK = reshape(KE(:)*((E0-Emin)*x(:)'.^penal+Emin),64*nelx*nely,1);
    K = sparse(iK, jK, sK); K = (K+K')/2;

    K(din, din) = K(din,din) +0.1;
    K(dout, dout) = K(dout,dout) +0.1;
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
    U(fixeddofs,:)= 0;
    U1 = U(:, 1);
    U2 = U(:, 2);
    c = full(U(dout,1));

    ce = reshape(sum((U1(edofMat) * KE) .* U2(edofMat), 2), nely, nelx);
    dc = penal*(E0-Emin)*x.^(penal-1).*ce;
  
% FILTERING OF SENSITIVITIES
    [dc]   = check(nelx,nely,rmin,x,dc);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    optTimer1 = tic;
    [x]    = OC(nelx,nely,x,volfrac,dc);    
    
    optTime = optTime + toc(optTimer1);
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% % PLOT DENSITIES  
%  colormap(gray); imagesc(1-x); axis equal; axis tight; axis off;pause(1e-6);  
end
totalTime = toc(t1);
disp(['total time: ', num2str(totalTime), 's']);
disp(['optimization time: ', num2str(optTime), 's']);
save('xm', 'x');
close all;
figure;
set(gcf, 'position', [200, 200, 400, 200])
colormap(gray); imagesc(1-x); axis equal; axis tight; axis off;pause(1e-6);

Emin = 1e-9;
level = 0.5;
levelHigh = 1;
levelLow = Emin;
xPhys = x;
while (true)
    x = xPhys;
    x(xPhys>=level) = 1;
    x(xPhys<level) = Emin;
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
[~, Uout]=FE(nelx,nely,xPhys,penal);
disp(Uout);
figure;
set(gcf, 'position', [200, 200, 400, 200])
colormap(gray); imagesc(1-xPhys); axis equal; axis tight; axis off;pause(1e-6);
end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.1;
while (l2-l1)/(l2+l1) > 1e-4 & l2 >1e-40;
  lmid = 0.5*(l2+l1);
  xnew = max(1e-3,max(x-move,min(1.,min(x+move,x.*(max(1e-10,-dc./lmid)).^0.3))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
end


%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U, Uout]=FE(nelx,nely,x,penal)
[KE] = lk; 
E0 = 1;
Emin = 1e-3;
F = sparse(2*(nely+1)*(nelx+1),2); U = sparse(2*(nely+1)*(nelx+1),2);
% K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
% for elx = 1:nelx
%   for ely = 1:nely
%     n1 = (nely+1)*(elx-1)+ely; 
%     n2 = (nely+1)* elx   +ely;
%     edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
%     K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
%   end
% end
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([-2 -1 2*nely+[0 1 2 3] 0 1],nelx*nely,1);
iK = reshape(kron(edofMat, ones(8, 1))', 64 * nelx * nely, 1);
jK = reshape(kron(edofMat, ones(1, 8))', 64 * nelx * nely, 1);
sK = reshape(KE(:)*((E0-Emin)*x(:)'.^penal+Emin),64*nelx*nely,1);
K = sparse(iK, jK, sK); K = (K+K')/2;
% DEFINE LOADS AND SUPPORTS (HALF FORCE INVERTER)
din = 1;
dout = 2*nelx*(nely+1)+1;
F(din, 1) = 1;
F(dout, 2) = -1;
K(din, din) = K(din,din) +0.1;
K(dout, dout) = K(dout,dout) +0.1;
fixeddofs   = union([2:2*(nely+1):2*(nelx+1)*(nely+1)],[2*(nely+1):-1:2*(nely+1)-5]);
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
Uout =full(U(dout,1));
end
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 1.; 
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
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
