% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 2 - Caso de estudio 2 - 
%   Inciso 3  
% Calcular un controlador que haga evolucionar al péndulo en el equilibrio estable, 
% partiendo de una condición inicial nula en el desplazamiento y 
% el ángulo en pi que termine en 2 metros evitando las oscilaciones de la masa m, 
% considerando que es una grúa. Una vez que delta=2 modificar a m a un 
% valor 10 veces mayor y volver al origen evitando oscilaciones.
% 
%   Inciso 4 
% Incorporar un observador para el caso en que sólo puedan medirse el desplazamiento
% delta y el ángulo phi, repetir las simulaciones para las condiciones anteriores 
% y graficar los resultados en gráficas superpuestas para el equilibrio estable.
%%
close all; clc; clear all;
m=.1;Fricc=0.1; long=1.6;g=9.8;M=1.5;
h=0.001;tiempo=(100/h);p_pp=0;tita_pp=0;
%Condiciones iniciales
alfa(1)=pi; color='m';

omega(1)=0; p_p(1)=0; u(1)=0; p(1)=0; i=1;
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0];
Mat_B=[0; 1/M; 0; 1/(long*M)];
Mat_C=[1 0 0 0]; %La salida es posición

% Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B ];%Matriz Controlabilidad
% rank(Mat_M);% verifico el rango de la matriz de controlabilidad

% %METODO LQR ------Controlador para m=.1 de 0 a 2 metros
Q = diag([.9 .2 .6 .2]);% Posicion Carro, Velocidad Carro, Angulo, Velocidad angular
R = 1000;
[K_m1,S,P] = lqr(Mat_A,Mat_B,Q,R);
G_m1=-inv(Mat_C*(eye(4)-inv(Mat_A+Mat_B*K_m1))*Mat_B);

% ------Controlador para m=1 de 2 a 0 metros
% Necesito agregar un integrador
Mat_Aa = [Mat_A zeros(4,1);-Mat_C 0];
Mat_Ba = [Mat_B; 0];
Mat_Ca=[Mat_C 0];

%  Controlador LQR
Q1 = diag([0.59 .01 .5 .006 .0001]);% Posicion Carro, Velocidad Carro, Angulo, Velocidad angular
R1 = 1000;
[K_m10,S_m10,P_m10] = lqr(Mat_Aa,Mat_Ba,Q1,R1);
G_m10= K_m10(1);


%OBSERVADOR
%redefino la matriz C para el observador
Mat_Co=[1 0 0 0;0 0 1 0]; %La salida es posición y ángulo

% % ------Controlador para m=.1 de 0 a 2 metros
Mat_A_O=Mat_A';
Mat_B_O=Mat_Co';
% Ko_m1=place(Mat_A_O,Mat_B_O,[-200 -300 -400 -500]);
% Qo = diag([1 .05 1 .05]); 
Qo=Mat_Co'*Mat_Co;
Ro = 100;
Ko_m1 = lqr(Mat_A_O,Mat_B_O,Qo,Ro);
% Go=Ko_m1(1); 
Ko=Ko_m1';
% Mat_M_Dual=[Mat_B_O Mat_A_O*Mat_B_O Mat_A_O^2*Mat_B_O Mat_A_O^3*Mat_B_O];%Matriz Controlabilidad
% eig(Mat_A_O'-Ko*Mat_C) %Verifico que todos los polos estén en el semiplano izquierdo

% % ------Controlador para m=1 de 2 a 0 metros
Mat_A_O=Mat_A';
Mat_B_O=Mat_Co';
Qo_1 = diag([1 .1 1 .1]); 
Ro_1 = [100 0;0 100];
[Ko_m10,So_m10,Po_m10] = lqr(Mat_A_O,Mat_B_O,Qo_1,Ro_1);
% Go=Ko_m10(1); 
% Ko=Ko_m10';
% Mat_M_Dual=[Mat_B_O Mat_A_O*Mat_B_O Mat_A_O^2*Mat_B_O Mat_A_O^3*Mat_B_O];%Matriz Controlabilidad
% eig(Mat_A_O'-Ko*Mat_C) %Verifico que todos los polos estén en el semiplano izquierdo

 
x_hat=[0 0;0 0;0 0]; %Inicializo el Observador

ref=2; 
psi(1)=0;
K=K_m1; G=G_m1;
while(i<(tiempo+1))
estado = [p(i); p_p(i); alfa(i); omega(i)]; % Posicion Carro, Velocidad Carro, Angulo, Velocidad angular
estado_obs = [p(i); p_p(i); alfa(i); omega(i)];
psi_p = ref - Mat_C* estado;
psi (i+1) = psi(i)+psi_p*h;


    if (i>=40000)
        ref=0;
        K = K_m10(1:4);
        Ki = -K_m10(5);
        G = G_m10;
        Go=Ko_m10(1); 
        Ko=Ko_m10';
%         u(i)=-K*estado+G*ref+Ki*psi(i+1); color='r' %LQR con integrador , sin observador
        u(i)=-K*x_hat+G*ref; color='b';% con observador
    else
%         u(i)=-K*estado+G*ref;  %color='r';% LQR sin integrador ,sin observador
        u(i)=-K*x_hat+G*ref; %color='b';% con observador
%     
    end
%  
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
subplot(3,2,1);plot(t,omega,color,'LineWidth',1.5);grid on; title('Velocidad ángulo');hold on;
% legend('LQR para referencia 2m y m=0.1');legend('boxoff');
legend('Sin Observador','Con Observador');legend('boxoff');
subplot(3,2,2);plot(t,alfa,color,'LineWidth',1.5);grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,p,color,'LineWidth',1.5);grid on;title('Posición grúa');hold on;
subplot(3,2,4);plot(t,p_p,color,'LineWidth',1.5);grid on;title('Velocidad grúa');hold on;
subplot(3,1,3);plot(t(1:end-1),u,color,'LineWidth',1.5);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;
% figure(2);hold on;
% subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
% subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;legend('Sin Observador','Con Observador');legend('boxoff');
