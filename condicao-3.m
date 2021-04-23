% Definição da Taxa de Decaimento
gama=0
% Declaração das Matrizes
A=[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
B=[0 0; 1 0; 0 0; 0 1];
C=[1 1 0 0; 0 0 1 1];
% Declaração da Matriz A considerando a Taxa de Decaimento
Ag=A+gama*eye(4)
% Declaração das Matrizes para o Complemento de Schur
I=eye(2);
eps=10^-4;
e2I=eps^2*I;
% Declaração das LMIs
setlmis([]);
P=lmivar(1,[4 1]);
Y=lmivar(2,[ 2 2]);
Z=lmivar(2,[2 2]);
% AgP+PAg’-BZC-C’Z’B’<0
lmiterm([1 1 1 P],Ag,1,’s’); % LMI #1: Ag*P+P*Ag’
lmiterm([1 1 1 Z],B,-C,’s’); % LMI #1: -B*Z*C-C’*Z’*B’
% P>0
lmiterm([-2 1 1 P],1,1); % LMI #2: P
% [e2I CP-YC; PC’-C’Y’ 1]>0
lmiterm([-3 1 1 0],e2I); % LMI #3: e2I
lmiterm([-3 2 1 P],1,C’); % LMI #3: P*C’
lmiterm([-3 2 1 -Y],C’,-1); % LMI #3: -C’*Y’
lmiterm([-3 2 2 0],1); % LMI #3: 1
manip=getlmis;
[tmin,xo]=feasp(manip);
% Resultado da Matriz P
Ps=dec2mat(manip,xo,P)
% Resultado da Matriz Y
Ys=dec2mat(manip,xo,Y)
% Resultado da Matriz Z
Zs=dec2mat(manip,xo,Z)
% Resultado da Matriz K
K=Zs*inv(Ys)
% Resultado dos polos do sistema para análise da estabilidade
eig(A-B*K*C)
