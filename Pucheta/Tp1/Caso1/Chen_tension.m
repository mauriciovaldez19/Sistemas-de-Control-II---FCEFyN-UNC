% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Ing. Pucheta
% Alumno: Valdez Benavidez, Mauricio Luciano
% Consigna 1 - Tp N° 1

%% Lectura de los datos de planilla excel
clear all; close all; clc;

valores=xlsread('Curvas_Medidas_RLC.xls'); 
t=valores(1:end,1);
I=valores(1:end,2);
Vc=valores(1:end,3);

%Grafico de tensión y corriente de los valores importados
figure(1)
plot(t,Vc, 'b' );title('Tension , V_t'); grid on;hold on; 

%Defino la entrada para la simulación posterior
t=linspace(0,0.1,1000);
u=linspace(0,0,1000);
vin=12;
ii=0;

for i=1:1000-1
 ii=ii+1;
 
 if ii<100
      u(i)=0;
 elseif ii>=100 && ii<=500
      u(i)=vin;
 else
      u(i)=vin*-1;
  end
end
figure(3)
plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on

% Aplicando el Metodo de Chen

%Se eligen los tres puntos según el artículo de Chen
t1_v=valores(126,1); y1_v=valores(126,3); 
t2_v=valores(151,1); y2_v=valores(151,3);
t3_v=valores(176,1); y3_v=valores(176,3);

figure(1)
plot([t1_v t2_v t3_v],[y1_v,y2_v,y3_v],'+');hold on;

%Ganancia seleccionada desde el grafico

K_v=12;     

%Defino las 3 k correspondientes a las 3 ecuaciones para los puntos tomados
k1_v= (y1_v/K_v) - 1;   
k2_v= (y2_v/K_v) - 1; 
k3_v= (y3_v/K_v) - 1;  


%Despejo de las ecuaciones alfa 1, alfa 2 y beta
b_v=4*(k1_v^3)*k3_v-3*(k1_v^2)*(k2_v^2)-4*(k2_v^3)+(k3_v^2)+6*k1_v*k2_v*k3_v; %EC 23 Chen
alfa1_v=(k1_v*k2_v+k3_v-sqrt(b_v))/(2*(k1_v^2+k2_v)); %EC 21 Chen
alfa2_v=(k1_v*k2_v+k3_v+sqrt(b_v))/(2*(k1_v^2+k2_v)); %EC 22 Chen
beta_v=(k1_v+alfa2_v)/(alfa1_v-alfa2_v);
%beta_v=(2*k1_v^3+3*k1_v*k2_v+k3_v-sqrt(b_v))/(sqrt(b_v)); %EC 24 Chen

%Sustituyendo EC 21 y EC 24 en EC 19 obtengo el cero y ambos polos
T1_v= -(t1_v-0.01)/log(alfa1_v);           
T2_v= -(t1_v-0.01)/log(alfa2_v);           
T3_v=(beta_v*(T1_v-T2_v))+T1_v;     

%Hago la nueva funcion de transferencia de tensión
G_v2=tf(K_v*[T3_v 1],conv([T1_v 1],[T2_v 1]));
G_v=tf(K_v,conv([T1_v 1],[T2_v 1]));

%Se compara la Tensión y la corriente

[y_G_v,t_G_v]=lsim(G_v,u/12,t);
figure(2)
plot(t_G_v,y_G_v,'k--'); title('Comparación de Tensión');
legend({'Vc(t) aproximada','Puntos seleccionados','Vc(t)real'},'Location','southeast')

%De la función de transferencia sacamos
% C=G_i.num{1}(2)
% L=(G_i.den{1}(1))/C
% R=(G_i.den{1}(2))/C


% L=0.125;
% R=223.14;
% Cap=2.1856e-5;
%Matrices
% A=[-R/L -1/L; 1/C 0];
% B=[1/L; 0];
% C=[1; 0]';
% D=0;
% %Definicion de la ecuación de estado y de salida (salida de corriente)
% G1=ss(A,B,C,D);
% [yout,yt]=lsim(G1,u,t);
% figure(3)
% plot(yt,yout,'b');grid on; 
% hold on;
% plot(valores(:,1),valores(:,2),'r');
% title('Comparación de corriente');
% legend({'i(t) aproximada','i(t) real'},'Location','southeast')
