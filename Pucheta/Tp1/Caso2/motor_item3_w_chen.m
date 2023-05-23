% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 1 - Caso de estudio 2 - 
%   Inciso 3 Parte 1.1
%   A partir de las curvas de mediciones de las variables graficadas en la Fig. 1-3, se requiere 
%   obtener el modelo del sistema considerando como entrada un escalón de 12V, como salida 
%   a la velocidad angular, y a partir de 0,1segundo se aplica un TL aproximado de 7,5 10-2
%   Nm. En el archivo Curvas_Medidas_Motor.xls están las mediciones, en la primer hoja
%   los valores y en la segunda los nombres. Se requiere obtener el modelo dinámico, para 
%   establecer las constantes de la corriente.

%% Lectura de los datos de planilla excel
clear all; close all; clc;

valores=xlsread('Curvas_Medidas_Motor_2023.xlsx'); 
tt=valores(1:end,1);
W=valores(1:end,2);
Ia=valores(1:end,3);
Vin=valores(1:end,4);
TL_=valores(1:end,5);
%Grafico de tensión y corriente de los valores importados
figure(1)
subplot(2,1,1);hold on
plot(tt,W, 'b' );title('Velocidad angular , \omega[rad/seg]'); grid on;hold on; 
subplot(2,1,2)
plot(tt,Vin, 'b' );title('Tensión,[V]');grid on;hold on;

%Defino la entrada para la simulación posterior
t_etapa=1e-7;
tF=.6;
u=linspace(0,0,(100e3-1));
ii=0;
for t=0:t_etapa:tF
ii=ii+1;
    u(ii)=0;
    if(t>=0.025)
        u(ii)=12;
    end
    
    if(t>=0.1501)       % Hago variar el valor de la tension de entrada
         u(ii)=-12;     % para simular el comporamiento de la tension con la Tl agregada
     end
end
t=0:t_etapa:tF;
% figure(4)
% plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on% 

 % Aplicando el Metodo de Chen
 
 %Se eligen los tres puntos según el artículo de Chen
t1=valores(1130,1); y1=(valores(1130,2)); 
t2=valores(2160,1); y2=(valores(2160,2));
t3=valores(3190,1); y3=(valores(3190,2));

figure(1)
subplot(2,1,1);hold on
plot([t1 t2 t3],[y1,y2,y3],'+');hold on;

%Ganancia seleccionada desde el grafico

K=(valores(15306,2))/12;
y1=y1/12;
y2=y2/12;
y3=y3/12;

%Defino las 3 k correspondientes a las 3 ecuaciones para los puntos tomados
k1= (y1/K) - 1;   
k2= (y2/K) - 1;   
k3= (y3/K) - 1;   


%Despejo de las ecuaciones alfa 1, alfa 2 y beta
b=4*(k1^3)*k3-3*(k1^2)*(k2^2)-4*(k2^3)+(k3^2)+6*k1*k2*k3; %EC 23 Chen
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2)); %EC 21 Chen
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2)); %EC 22 Chen
beta=(k1+alfa2)/(alfa1-alfa2); %EC 20 Chen
%beta=(k2+alfa2^2)/(alfa1^2-alfa2^2); %EC 20 Chen
%beta=(k3+alfa2^3)/(alfa1^3-alfa2^3); %EC 20 Chen

%beta=(2*(k1^3)+3*k1*k2+k3-sqrt(b))/(sqrt(b)); %EC 24 Chen

%Sustituyendo EC 21 y EC 24 en EC 19 obtengo el cero y ambos polos
T1= -(t1-0.025)/log(alfa1);      
T2= -(t1-0.025)/log(alfa2);      
T3=(beta*(T1-T2))+T1;            


%Hago la nueva funcion de transferencia de tensión
G_w=tf(K*[T3 1],conv([T1 1],[T2 1]));

[y_G_w,t_G_w]=lsim(G_w,u,t);

figure(2)
plot(tt,W, 'k' ); grid on; hold on;
plot(t_G_w,y_G_w,'r'); title('G_w (s) obtenida con método de Chen vs \omega de tabla');
legend({'w_t de excel','G_w (s) obtenida con método de Chen'},'Location','southeast')


