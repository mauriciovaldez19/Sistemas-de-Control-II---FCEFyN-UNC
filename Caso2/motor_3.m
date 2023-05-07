% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 1 - Caso de estudio 2 - 
%   Inciso 3 Parte 2
%   A partir de las curvas de mediciones de las variables graficadas en la Fig. 1-3, se requiere 
%   obtener el modelo del sistema considerando como entrada un escalón de 12V, como salida 
%   a la velocidad angular, y a partir de 0,1segundo se aplica un TL aproximado de 7,5 10-2
%   Nm. En el archivo Curvas_Medidas_Motor.xls están las mediciones, en la primer hoja
%   los valores y en la segunda los nombres. Se requiere obtener el modelo dinámico, para 
%   establecer las constantes de la corriente.
%%
clc;clear;close all;
X=-[0; 0; 0; 0];ii=0;t_etapa=1e-7;tF=.6;
color_='r';
Ts=t_etapa;
%u=12;
for t=0:t_etapa:tF
ii=ii+1;k=ii+2;
Tl=0; u=0;
    if(t>=0.025)
        u=12;
    end
%     if(t>=0.125)        %tiempo indicado por consigna
%         Tl=7.5e-2;    %Tl indicado por consigna
%    end
     if(t>=0.1501)  % tiempo tabla excel
%          Tl=7.5e-2; % Tl indicado por la consigna
         Tl=4.5e-2; % Tl que se ajusta a la rta de los datos del excel
     end
     
X=modmotorpunto2(t_etapa, X, u,Tl);
x1(ii)=X(1);%Omega
x2(ii)=X(2);%wp
x3(ii)=X(3);%ia
x4(ii)=X(4);%tita
acc(ii)=u;
end
t=0:t_etapa:tF;

%Cargo valores para la Comparacion
valores=xlsread('Curvas_Medidas_Motor_2023.xlsx'); 
tt=valores(1:end,1);
W=valores(1:end,2);


figure(1)
subplot(2,1,1);hold on;
plot(tt,W, 'g' );title('Velocidad angular , W[rad/seg]'); grid on;hold on; 
plot(t,x1,color_);title('Salida y, \omega_t');hold on;
legend({'w de excel','w aproximada'},'Location','southeast')

subplot(2,1,2);hold on;
plot(t,acc,'r');
xlabel('Tiempo [Seg.]');hold on;




%motor
function [X]=modmotorpunto2(t_etapa, xant, accion,Tl)
Laa=506.94e-6; J=3.963e-6;Ra=99.7224;B=0;Ki=16.52;Km=(1/16.52);
Va=accion;
h=1e-7;
TL=Tl;   
omega= xant(1);     
wp= xant(2);        
ia=xant(3);
tita = xant(4);
for ii=1:t_etapa/h
wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
wp=wp+h*wpp;
wp=wp-(TL/J);
ia=0;
omega = omega + h*wp;
tita = tita + h*omega;
end

X=[omega,wp,ia,tita];
end