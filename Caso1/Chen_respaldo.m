% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Ing. Pucheta
% Alumno: Valdez Benavidez, Mauricio Luciano
% Consigna 1 - Tp N° 1

clear all; close all; clc;

valores=xlsread('Curvas_Medidas_RLC.xls'); 
t=valores(1:end,1);
I=valores(1:end,2);
Vc=valores(1:end,3);

%Grafico de tensión
figure(1)
plot(t,Vc,'m');grid on; title('Tensión');hold on;

%Se eligen los tres puntos según el artículo de Chen
t1=valores(106,1); y1=valores(106,3); 
t2=valores(116,1); y2=valores(116,3);
t3=valores(126,1); y3=valores(126,3);
plot([t1 t2 t3],[y1,y2,y3],'ok')

%Definición de la entrada
u=zeros(1,1000);
paso=0.1/1000;
t=0:paso:(0.1-paso);
signo=true;
for i=100:1:1000
 if mod(i,500)==0
 signo=not(signo);
 end
 if signo==1
 u(1,i)=12;
 end
 if signo==0
 u(1,i)=-12;
 end
end
%Ganancia de corriente
k=12; 
%Defino las 3 k correspondientes a las 3 ecuaciones para los puntos tomados
k1= (y1/k) - 1;
k2= (y2/k) - 1;
k3= (y3/k) - 1;
%Despejo de las ecuaciones alfa 1, alfa 2 y beta
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3; %EC 23 Chen
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2)); %EC 21 Chen
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2)); %EC 22 Chen
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b)); %EC 20 Chen
%Sustituyendo EC 21 y EC 24 en EC 19 obtengo el cero y ambos polos
T1= -0.001/log(alfa1);
T2= -0.001/log(alfa2);
T3= (beta*(T1-T2))+T1; %No hay cero, no se usa
%Hago la nueva funcion de transferencia de tensión
G=tf(k,conv([T1 1],[T2 1]))
%Se compara la Tensión
[yaprox,taprox]=lsim(G,u/12,t);
hold on;
plot(taprox,yaprox,'k--'); title('Comparación de Tensión');
legend({'Vc(t) aproximada','Puntos seleccionados','Vc(t)real'},'Location','southeast')
%De la función de transferencia sacamos
L=0.125;
R=223.14;
Cap=2.1856e-5;
%Matrices
A=[-R/L -1/L; 1/Cap 0];
B=[1/L; 0];
C=[1; 0]';
D=0;
%Definicion de la ecuación de estado y de salida (salida de corriente)
G1=ss(A,B,C,D);
[yout,yt]=lsim(G1,u,t);
figure
plot(yt,yout,'b');grid on; 
hold on;
plot(valores(:,1),valores(:,2),'r');
title('Comparación de corriente');
legend({'i(t) aproximada','i(t) real'},'Location','southeast')
