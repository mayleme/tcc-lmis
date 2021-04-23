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
N=lmivar(2,[2 2]);
M=lmivar(2,[2 2]);
% P*Ag+Ag’*P-B*N*C-C’*N’*B’<0
  lmiterm([1 1 1 P],1,Ag,’s’); % LMI #1: P*Ag+Ag’*P
lmiterm([1 1 1 N],B,-C,’s’); % LMI #1: -B*N*C-C’*N’*B’
% P>0
lmiterm([-2 1 1 P],1,1); % LMI #2: P
% [e2I B’*P-M*B’; P*B-B*M’ 1]>0
lmiterm([-3 1 1 0],e2I); % LMI #3: e2I
lmiterm([-3 2 1 P],1,B); % LMI #3: P*B
lmiterm([-3 2 1 0],-B*M’); % LMI #3: -B*M’
lmiterm([-3 2 2 0],1); % LMI #3: 1
manip=getlmis;
[tmin,xo]=feasp(manip);
% Resultado da Matriz P
Ps=dec2mat(manip,xo,P)
% Resultado da Matriz N
Ns=dec2mat(manip,xo,N)
% Resultado da Matriz M
Ms=dec2mat(manip,xo,M)
% Resultado da Matriz K
K=inv(Ms’)*Ns
% Resultado dos polos do sistema para análise da estabilidade
eig(A-B*K*C)
