% Definição da Taxa de Decaimento
gama=0
% Declaração das Matrizes
A=[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
B=[0 0; 1 0; 0 0; 0 1];
C=[1 1 0 0; 0 0 1 1];
In=eye(4);
d=1
dIn=d*In;
% Declaração da Matriz A considerando a Taxa de Decaimento
Ag=A+gama*eye(4)
% Declaração das Matrizes para o Complemento de Schur
eps=10^-4;
e2I=eps^2*In;
% Declaração das LMIs
setlmis([]);
P=lmivar(1,[4 1]);
K=lmivar(1,[2 1]);
Fn=lmivar(2,[2 2]);
%[-2P A’-C’K’B’+P+dIn; A-BKC+P+dIn -2dIn]<0
lmiterm([1 1 1 P],.5*2,-1,’s’); % LMI #1: -2*P
lmiterm([1 2 1 K],B,-C); % LMI #1: -B*K*C
lmiterm([1 2 1 P],1,1); % LMI #1: P
lmiterm([1 2 1 0],Ag+dIn); % LMI #1: Ag+dIn
lmiterm([1 2 2 0],-2*dIn); % LMI #1: -2*dIn
%P>0
lmiterm([-2 1 1 P],1,1); % LMI #2: P
% [e2I CP-YC; PC’-C’Y’ 1]>0
lmiterm([-3 1 1 0],e2I); % LMI #3: e2I
lmiterm([-3 2 1 P],B’,1); % LMI #3: B’*P
lmiterm([-3 2 1 Fn],1,-C); % LMI #3: -Fn*C
lmiterm([-3 2 2 0],1); % LMI #3: 1
manip=getlmis;
[tmin,xo]=feasp(manip);
% Resultado da Matriz P
Ps=dec2mat(manip,xo,P)
% Resultado da Matriz K
Ks=dec2mat(manip,xo,K)
% Resultado da Matriz Fn
Fs=dec2mat(manip,xo,Fn)
%Resultado da Matriz F
F=(1/d)*Fn
% Resultado dos polos do sistema para análise da estabilidade
eig(A-B*Ks*C)
