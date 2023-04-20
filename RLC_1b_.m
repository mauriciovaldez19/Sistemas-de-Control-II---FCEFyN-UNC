% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Ing. Pucheta
% Alumno: Valdez Benavidez, Mauricio Luciano
% Consigna 1 - Tp N° 1

clear all; close all; clc;
R=5.6e3; L=10e-6; C=100e-9; %Asignación de valores para item 2   

A=[-R/L -1/L; 1/C 0]; %matriz de estados
B=[1/L; 0]; %matriz de entrada
C=[R 0]; %matriz de salida
D=[0]; %matriz de transmision directa

[num,den]=ss2tf(A,B,C,D); %obtengo la FT a partir del espacio de estados dado y luego los polos
G=tf(num,den); polos=roots(den);

tR=log(0.95)/polos(1); %intervalo de integracion -> dinámica mas rápida, tomo 1e-11 por cuestiones de simplicidad
tR=9e-10;
tL=log(0.05)/polos(2); %tiempo de simulacion para entrada invariable -> dinámica mas lenta, tomo 5-e3
tL=9e-3;

paso=tL/tR;
t=linspace(0,tL,paso);
I=zeros(1,paso); Vc=zeros(1,paso); u=linspace(0,0,paso);
vin=12;
%condiciones iniciales
I(1)=0; Vc(1)=0; u(1)=vin;

%punto de operacion
Il(1)=0;
Vcl(1)=0;
x=[I(1) Vc(1)]';
y(1)=0;
Xop=[0 0]';
ii=0;

for i=1:paso-1
 ii=ii+tR;
 if (ii>=1e-3)
 ii=0;
 vin=vin*-1;
 end
 u(i)=vin;
 %Variables de sistema lineal
 xp=A*(x-Xop)+B*u(i);
 x=x+xp*tR;
 Y=C*x;
 y(i+1)=Y(1);
 Il(i+1)=x(1);
 Vcl(i+1)=x(2);
end
figure(1)
%Grafica I, Vc y u
subplot(3,1,1);
plot(t,Il, 'b' );title('Corriente , i_t'); grid on; 
subplot(3,1,2);
plot(t,Vcl, 'r' );title('Tensión Capacitor , Vc_t');grid on
subplot(3,1,3); 
plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on
