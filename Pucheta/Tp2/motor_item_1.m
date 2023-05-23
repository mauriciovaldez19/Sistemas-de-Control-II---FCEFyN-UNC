% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 2 - Caso de estudio 1 - 
%   Inciso 1  


%   
%%
clear all; close all; clc;
Laa=5e-3; J=0.004; Ra=0.2; B=0.005; Ki=6.5e-5;Km=0.055;

A=[-Ra/Laa -Km/Laa 0; Ki/J -B/J 0; 0 1 0];   %matriz de estados
B=[1/Laa 0; 0 -1/J; 0 0];             %matriz de entrada
C=[0 0 1];                %matriz de salida

polos=eig(A);

tR=log(0.95)/polos(3); %dinámica mas rápida
tR=tR/10;
tL=log(0.05)/polos(2); %dinámica mas lenta
tL=tL*5;

paso=tL/tR;
t=linspace(0,tL,paso);
u=linspace(0,0,paso);
vin=12;


%Condiciones iniciales
Il(1)=0;
Vcl(1)=0;
x=[Il(1) Vcl(1)]';
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
%Grafico de Il, Vc y Vin
figure(1)
subplot(3,1,1);
plot(t,Il, 'b' );title('Corriente , i_t'); grid on; 
subplot(3,1,2);
plot(t,Vcl, 'r' );title('Tensión Capacitor , Vc_t');grid on
subplot(3,1,3); 
plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on

%Grafico de Il y Vin
figure(2)
subplot(2,1,1);
plot(t,Il, 'b' );title('Corriente , i_t'); grid on; 
subplot(2,1,2); 
plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on

%Grafico de Vc y Vin
figure(3)
subplot(2,1,1);
plot(t,Vcl, 'r' );title('Tensión Capacitor , Vc_t');grid on
subplot(2,1,2); 
plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on

%Grafico de Il y Vc
figure(4)
subplot(2,1,1);
plot(t,Il, 'b' );title('Corriente , i_t'); grid on; 
subplot(2,1,2); 
plot(t,Vcl, 'r' );title('Tensión Capacitor , Vc_t');grid on
