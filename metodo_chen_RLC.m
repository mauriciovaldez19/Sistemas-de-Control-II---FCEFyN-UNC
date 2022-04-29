clear all;close all;clc
syms T1 T2 T3 k1 k2 k3 alfa1 alfa2 beta T1_I T2_I T3_I k1_I k2_I k3_I alfa1_I alfa2_I beta_I
% Importo las mediciones dadas
data=xlsread('Curvas_Medidas_RLC.xls');

t = data(1:end,1);
I = data(1:end,2);
Vc = data(1:end,3);

pos= 110; 

%% Calculo FT para la tension del capacitor
% Tomo 3 puntos (t1,y(t1),(2t1,y(2t1),(3t1,y(3t1), defino la variable pos, p/ selecc distintas posiciones de la tabla y realizar pruebas
t1 = t(pos); y1 = Vc(pos);
t2 = t(2*pos); y2 = Vc(2*pos);
t3 = t(3*pos); y3 = Vc(3*pos);

% Valor de K en regimen, elijo el último de la tabla
K=abs(Vc(end));
% Defino k1, k2 y k3 correspondientes a las 3 ecuaciones para los 3 puntos tomados
k1 = (y1/K)-1;
k2 = (y2/K)-1;
k3 = (y3/K)-1;
%Despejo de las ecuaciones alfa 1, alfa 2 y beta

eqn1=alfa1*alfa2==((k2^2)-(k1*k3))/((k1^2)+k2);
eqn2=alfa1+alfa2==(k3+(k1*k2))/((k1^2)+k2);
eqns = [eqn1 eqn2];
s=solve(eqns,[alfa1 alfa2],'ReturnConditions',true);
simplify(s.alfa1);
simplify(s.alfa2);
alfa1=s.alfa1(1);
alfa2=s.alfa2(2);

be=(4*k1^3*k3)-(3*k1^2*k2^2)-(4*k2^3)+(k3^2)+(6*k1*k2*k3); 
alfa1=((k1*k2)+k3-sqrt(be))/(2*(k1^2+k2)); 
alfa2=((k1*k2)+k3+sqrt(be))/(2*(k1^2+k2)); 

beta = (k1+alfa2)/(alfa1-alfa2); 

T1_ang=-t1/log(alfa1);
T2_ang=-t1/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;

% Escalón unitario para la simulación
escalon = stepDataOptions;
escalon.StepAmplitude,1;
% Defino la funcion de transferencia y la dibujo
G_Vc=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))

%% Calculo FT para la Corriente del inductor
% Tomo 3 puntos (t1,y(t1),(2t1,y(2t1),(3t1,y(3t1), defino la variable pos, p/ selecc distintas posiciones de la tabla y realizar pruebas

t1_I = t(pos); y1_I = I(pos);
t2_I = t(2*pos); y2_I = I(2*pos);
t3_I = t(3*pos); y3_I = I(3*pos);

% Valor de K en regimen, elijo el último de la tabla
K_I=abs(I(end));
% Defino k1, k2 y k3 correspondientes a las 3 ecuaciones para los 3 puntos tomados
k1_I = (y1_I/K_I)-1;
k2_I = (y2_I/K_I)-1;
k3_I = (y3_I/K_I)-1;
%Despejo de las ecuaciones alfa 1, alfa 2 y beta

eqn1_I=alfa1_I*alfa2_I==((k2_I^2)-(k1_I*k3_I))/((k1_I^2)+k2_I);
eqn2_I=alfa1_I+alfa2_I==(k3_I+(k1_I*k2_I))/((k1_I^2)+k2_I);
eqns_I = [eqn1_I eqn2_I];
s_I=solve(eqns_I,[alfa1_I alfa2_I],'ReturnConditions',true);
simplify(s_I.alfa1_I);
simplify(s_I.alfa2_I);
alfa1_I=s_I.alfa1_I(1);
alfa2_I=s_I.alfa2_I(2);

be_I=(4*k1_I^3*k3_I)-(3*k1_I^2*k2_I^2)-(4*k2_I^3)+(k3_I^2)+(6*k1_I*k2_I*k3_I); 
alfa1_I=((k1_I*k2_I)+k3_I-sqrt(be_I))/(2*(k1_I^2+k2_I)); 
alfa2_I=((k1_I*k2_I)+k3_I+sqrt(be_I))/(2*(k1_I^2+k2_I)); 

beta_I = (k1_I+alfa2_I)/(alfa1_I-alfa2_I); 

T1_ang_I=-t1_I/log(alfa1_I);
T2_ang_I=-t1_I/log(alfa2_I);
T3_ang_I=beta_I*(T1_ang_I-T2_ang_I)+T1_ang_I;

% Escalón unitario para la simulación
escalon = stepDataOptions;
escalon.StepAmplitude,1;
% Defino la funcion de transferencia y la dibujo
G_I=tf(K_I*[T3_ang_I 1],conv([T1_ang_I 1],[T2_ang_I 1]))

%% Despejo los valores de RLC
C=G_I.num{1}(2)
L=(G_I.den{1}(1))/C
R=(G_I.den{1}(2))/C

%Calculo la nueva G con los valores de RLC con la ganancia K de la G_Vc
G_aprox=tf(K,[L*C R*C 1]);


%% Ploteo
subplot(5,1,1);plot(t,Vc);grid on; title('Tensón en Capacitor');
hold on;plot([t1 t2 t3],[y1,y2,y3],'+')
subplot(5,1,2);plot(t,I);grid on; title('Corriente');
hold on;plot([t1_I t2_I t3_I],[y1_I,y2_I,y3_I],'O')

subplot(5,1,3);step(G_Vc,escalon,'r'),hold on;grid on;title('Tensión Estimada');
subplot(5,1,4);step(G_I,escalon,'y'),hold on;grid on;title('Corriente Estimada');

subplot(5,1,2);step(G_aprox,escalon,'g'),hold on;grid on;title('Tension Aproximada');
