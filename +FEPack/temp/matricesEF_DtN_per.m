function [M,K] = matricesEF_DtN_per(const,func,maillage)

nDoF = maillage.nDoF;

II = zeros(maillage.nTri*9,1);
JJ = zeros(maillage.nTri*9,1);
VV_K = zeros(maillage.nTri*9,1);
VV_M = zeros(maillage.nTri*9,1);
index_II = [1 2 3 1 2 3 1 2 3];
index_JJ = [1 1 1 2 2 2 3 3 3];

coeff = func.per;
delta = const.delta;

for t = 1:maillage.nTri
    
    S1 = maillage.Coorneu(maillage.Numtri(t,1),:);
    S2 = maillage.Coorneu(maillage.Numtri(t,2),:);
    S3 = maillage.Coorneu(maillage.Numtri(t,3),:);
    
    % calcul des matrices elementaires du triangle t
    
    [Mloc,Kloc] = mat_elem_DtN_per(S1, S2, S3, coeff, delta);
    
    index = (9*(t-1)+1):(9*(t-1)+9);
    II(index) = maillage.Numtri(t,index_II);
    JJ(index) = maillage.Numtri(t,index_JJ);
    VV_K(index) = Kloc(:);
    VV_M(index) = Mloc(:);

end


M = sparse(II,JJ,VV_M,nDoF,nDoF);
K = sparse(II,JJ,VV_K,nDoF,nDoF);

