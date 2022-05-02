%%%%%% A 128 LINE TOPOLOGY OPTIMIZATION CODE BY Yuan Liang Gengdong Cheng %%%
function DVTOPCRA(nelx,nely,volfrac,penal,rmin,beta)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% Parameters
vrf = 0.98;        %Volume reduction factor
w1 = 1E-8;         %Error tolerance of Canonical Approximate method
w2 =5E-4;         %Threshold for reducting volume
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1), 2*(nely+1)*(nelx+1));
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = 1;
            end
        end
    end
end
H1 = sparse(iH,jH,sH);
HS1 = sum(H1,2);
N = nely*nelx;
%% INITIALIZE ITERATION
x_iteration = ones(nely,nelx);
Vr = vrf*sum(x_iteration(:));
Lambda = 1e-3;      %Initial multiplier for the material usage constraint
ID = ones(nely,nelx);
FEMnumber = 0;
Compliance = zeros(2,1);
VolumeFrac = zeros(2,1);
pareto = 0;
%% Structural analysis and sensitivity annalysis
sK = reshape(KE(:)*(Emin+x_iteration(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
if penal <= 1  dc = -x_iteration.*ce; else   dc = -penal*(x_iteration.^(penal-1)).*ce; end
dc(:) = H1*(dc(:)./HS1);
dc = dc./max(max(abs(dc)));
c_iteration = full(0.5*dot(F,U));
flag=0;
%% Start iterate
while 1         %Start reduce volume fraction
    x_old = x_iteration;
    %Implement Canonical Ralaxation algorithm
    while 1        
        s = Lambda.*ID+dc;
        index=abs(s) <= 1E-8;
        s(index) = 1E-8*sign(s(index));
        delta = beta^2/27;
        r = power(delta,-1/3).*power(2.*s.^2-delta.*ID+2.*sqrt(s.^2.*(s.^2-delta.*ID)),1/3);
        r_con = conj(r);
        sigma = (1/6)*beta.*(-ID+r+r_con);
        Lambda = (sum(sum(ID-dc./sigma))-2*Vr)/sum(sum(1./sigma));
        x_new = (1/2).*(ID-(Lambda.*ID+dc)./sigma);
        c_old = sum(sum(x_old.*dc));
        c_new = sum(sum(x_new.*dc));
        if abs(c_old-c_new)/abs(c_old) <= w1%&&sum(sum(x_new)) <= Vr
            break;   %Satisfy convergence criterion
        else
            x_old = x_new;    %Renew design variable
        end
    end
    x_new(x_new<=0) = 0.001;    
    x_new(x_new>=1) = 1;
    x_new=round(x_new);  %Eliminate numerical error
    figure(1); colormap(gray); imagesc(1-x_new); caxis([0 1]); axis equal; axis off; drawnow;
    sK = reshape(KE(:)*(Emin+x_new(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    if penal <=1 dc = -x_new.*ce; else  dc = -penal*(x_new.^(penal-1)).*ce;    end
    dc(:) = H1*(dc(:)./HS1);
    dc = dc./max(max(abs(dc)));
    c_new = full(0.5*dot(F,U));
    change = abs(c_new-c_iteration)/c_iteration;
    FEMnumber = FEMnumber+1;
    fprintf(' It.:%5i Obj.:%7.4f Vol.:%6.4f ch.:%8.4f\n',FEMnumber,2*c_new,mean(x_new(:)),change);
    if abs(c_new-c_iteration)/c_iteration > w2    %Process for current volume is not convergent
        x_iteration = x_new;
        c_iteration = c_new;
    else
        if flag==1  break;  end
        x_iteration = x_new;
        pareto = pareto+1;
        Compliance(pareto) = c_new;
        VolumeFrac(pareto) = mean(x_new(:));
        if Vr/N <= volfrac
            Vr =(max(volfrac,Vr/N))*N;
            flag=1;
            continue;
        else
            Vr = Vr*vrf;    %Reduce material usage
        end
    end
end
save('x', 'x_iteration');
figure(2); plot(VolumeFrac,Compliance,'b-o','LineWidth',2); ...
    xlabel('Volume fraction','Fontname', 'Times New Roman','FontSize',12);...
    ylabel('Structural Compliance','Fontname', 'Times New Roman','FontSize',12);...
    title('Pareto Frontier'); %Output Pareto Frontier
end