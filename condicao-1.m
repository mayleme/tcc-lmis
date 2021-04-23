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
R=lmivar(2,[2 2]);
F=lmivar(2,[2 2]);
% PAg+Ag’P-C’RC-C’R’C<0
lmiterm([1 1 1 P],1,Ag,’s’); % LMI #1: P*Ag+Ag’*P
lmiterm([1 1 1 R],C’,-C,’s’); % LMI #1: -C’*R*C-C’*R’*C
% P>0
lmiterm([-2 1 1 P],1,1); % LMI #2: P
% [e2I B’P-FC; PB-C’F’ 1]>0
lmiterm([-3 1 1 0],e2I); % LMI #3: e2I
lmiterm([-3 2 1 P],1,B); % LMI #3: P*B
lmiterm([-3 2 1 -F],C’,-1); % LMI #3: -C’*F’
lmiterm([-3 2 2 0],1); % LMI #3: 1
manip=getlmis;
[tmin,xo]=feasp(manip);
% Resultado da Matriz P
Ps=dec2mat(manip,xo,P)
% Resultado da Matriz R
Rs=dec2mat(manip,xo,R)
% Resultado da Matriz F
Fs=dec2mat(manip,xo,F)
% Resultado da Matriz K
K=inv(Fs’)*Rs
% Resultado dos polos do sistema para análise da estabilidade
eig(A-B*K*C)
