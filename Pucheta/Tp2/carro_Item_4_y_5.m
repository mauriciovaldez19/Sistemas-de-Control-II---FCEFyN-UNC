% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N� 2 - Caso de estudio 2 - 
%   Inciso 1  
% Calcular un controlador que haga evolucionar al p�ndulo en el equilibrio inestable, 
% partiendo de una condici�n inicial nula en el desplazamiento y termine en -10 metros manteniendo 
% la vertical. Determinar el �ngulo m�ximo que puede alejarse de la vertical en t=0 para que el sistema 
% cumpla el objetivo de control.
%   Inciso2
% Incorporar un observador para el caso en que s�lo puedan medirse el desplazamiento y el �ngulo, 
% repetir las simulaciones para las condiciones anteriores y graficar los resultados en 
% gr�ficas superpuestas.
%%
close all; clc; clear all;
m=.1;Fricc=0.1; long=1.6;g=9.8;M=1.5;
h=0.001;tiempo=(50/h);p_pp=0;tita_pp=0;
%Condiciones iniciales
alfa(1)=.1; color='r';
% alfa(1)=.4; color='b';
% alfa(1)=.8; color='m';
% alfa(1)=1.2; color='k'; %Angulo inicial m�ximo para LQR sin obs
% alfa(1)=1.3; color='g'; %Angulo inicial para inestabilizar el sistema para LQR sin obs
omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;indice=0;
%Versi�n linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_B=[0; 1/M; 0; -1/(long*M)];
Mat_C=[1 0 1 0]; %La salida es posici�n y �ngulo

% Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B ];%Matriz Controlabilidad
% rank(Mat_M);% verifico el rango de la matriz de controlabilidad

%METODO LQR
Q = diag([10 1 10 1]);% Posicion Carro, Velocidad Carro, Angulo, Velocidad angular
% Q=Mat_C'*Mat_C;
R = 1000;
[K,S,P] = lqr(Mat_A,Mat_B,Q,R);
G=-inv(Mat_C*(eye(4)-inv(Mat_A+Mat_B*K))*Mat_B); %%Es igual a Klqr(1)


%OBSERVADOR
Mat_A_O=Mat_A';
Mat_B_O=Mat_C';
Qo = diag([1000 50 1000 50]);
Ro = 10;
[Ko,So,Po] = lqr(Mat_A_O,Mat_B_O,Qo,Ro);
Go=Ko(1); %%Es igual a Klqr(1)
Ko=Ko';
% Mat_M_Dual=[Mat_B_O Mat_A_O*Mat_B_O Mat_A_O^2*Mat_B_O Mat_A_O^3*Mat_B_O];%Matriz Controlabilidad
% eig(Mat_A_O'-Ko*Mat_C) %Verifico que todos los polos est�n en el semiplano izquierdo
x_hat=[0;0;0;0]; %Inicializo el Observador

ref=-10; psi(1)=0;

while(i<(tiempo+1))
 estado=[p(i); p_p(i); alfa(i); omega(i)]; % Posicion Carro, Velocidad Carro, Angulo, Velocidad angular
 estado_obs=[p(i); p_p(i); alfa(i); omega(i)];

u(i)=-K*estado+G*ref; %color='r';% lqr sin observador
% u(i)=-K*x_hat+G*ref; color='m';% con observador

p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))-Fricc*p_p(i));
tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
 p_p(i+1)=p_p(i)+h*p_pp;
 p(i+1)=p(i)+h*p_p(i);
 omega(i+1)=omega(i)+h*tita_pp;
 alfa(i+1)=alfa(i)+h*omega(i);
 y_sal(i)=Mat_C*estado;
 %________OBSERVADOR__________
 y_sal_O(i)=Mat_C*x_hat;
 y_sal(i)=Mat_C*estado;
 x_hatp=Mat_A*x_hat+Mat_B*u(i)+Ko*(y_sal(i)-y_sal_O(i));
 x_hat=x_hat+h*x_hatp;
 i=i+1;
end
figure(1);hold on; t=1:i;t=t*h;
subplot(3,2,1);plot(t,omega,color,'LineWidth',1.5);grid on; title('Velocidad �ngulo');hold on;
% legend('\theta = 0.1','\theta = 0.4','\theta = 0.8','\theta = 1.2');legend('boxoff');
legend('Sin Observador','Con Observador');legend('boxoff');
subplot(3,2,2);plot(t,alfa,color,'LineWidth',1.5);grid on;title('�ngulo');hold on;
subplot(3,2,3); plot(t,p,color,'LineWidth',1.5);grid on;title('Posici�n carro');hold on;
subplot(3,2,4);plot(t,p_p,color,'LineWidth',1.5);grid on;title('Velocidad carro');hold on;
subplot(3,1,3);plot(t(1:end-1),u,color,'LineWidth',1.5);grid on;title('Acci�n de control');xlabel('Tiempo en Seg.');hold on;
% figure(2);hold on;
% subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('�ngulo');ylabel('Velocidad angular');hold on;
% subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;legend('Sin Observador','Con Observador');legend('boxoff');
