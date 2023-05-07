% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 1 - Caso de estudio 2 - 
%   Inciso 1  
%   Obtener el torque máximo que puede soportar el motor modelado mediante las Ecs. (5)
%   (6) y (7) cuando se lo alimenta con 12V, graficando para 5 segundos de tiempo la
%   velocidad angular y corriente ia.

%   Inciso 2
%   Mostrar simulaciones de 5 segundos que permitan observar la corriente ia en todo 
%   momento y establecer su valor máximo como para dimensionar dispositivos electrónicos.

%%
clc;clear;close all;
X=-[0; 0; 0; 0];ii=0;t_etapa=1e-7;tF=5;
color_='b';
Ts=t_etapa;
u=12;
for t=0:t_etapa:tF
ii=ii+1;k=ii+2;
X=modmotorpunto2(t_etapa, X, u);
x1(ii)=X(1);%Omega
x2(ii)=X(2);%wp
x3(ii)=X(3);%ia
x4(ii)=X(4);%tita
acc(ii)=u;
end
t=0:t_etapa:tF;
subplot(3,1,1);hold on;
plot(t,x1,color_);title('Salida y, \omega_t');
subplot(3,1,2);hold on;
plot(t,x3,'k');title('Corriente de salida, i_a');
subplot(3,1,3);hold on;
plot(t,acc,'b');title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');

%motor
function [X]=modmotorpunto2(t_etapa, xant, accion)
Laa=366e-6; J=5e-9;Ra=55.6;B=0;Ki=6.49e-3;Km=6.53e-3;
Va=accion;
h=1e-7;
TL=2.1278125e-5;    %se llego a ese valor luego de varias iteraciones.
omega= xant(1);     %Se partió desde 1 y se fue bajando dividiendo por 10
wp= xant(2);        %cada vez hasta que w fue + y luego busco algo de precisión
ia=xant(3);
tita = xant(4);
for ii=1:t_etapa/h
wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
iap=(-Ra*ia-Km*omega+Va)/Laa;
wp=wp+h*wpp;
wp=wp-((1/J)*TL);
ia=ia+iap*h;
omega = omega + h*wp;
tita = tita + h*omega;
end
X=[omega,wp,ia,tita];
end