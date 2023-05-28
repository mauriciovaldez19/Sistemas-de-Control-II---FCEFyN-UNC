% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Valdez Benavidez, Mauricio Luciano
% Tp N° 2 - Caso de estudio 2 - 
%   Inciso 1  
% Calcular un controlador que haga evolucionar al péndulo en el equilibrio inestable, 
% partiendo de una condición inicial nula en el desplazamiento y termine en -10 metros manteniendo 
% la vertical. Determinar el ángulo máximo que puede alejarse de la vertical en t=0 para que el sistema 
% cumpla el objetivo de control.
%%
clc;clear all;

%EQUILIBRIO INESTABLE

%PARAMETROS
m = 0.1; Fricc = 0.1; long = 1.6; g = 9.8; M = 1.5;

%MATRICES
Mat_A = [0 1 0 0;                               %X1 = delta
         0 -Fricc/M -m*g/M 0;                   %X2 = delta_p
         0 0 0 1;                               %X3 = phi
         0 Fricc/(long*M) g*(m+M)/(long*M) 0];  %X4 = phi_p
    
Mat_B = [0; 1/M; 0; -1/(long*M)];
Mat_C = [1 0 0 0] ;

% %     tiempo de integracion y tiempo de simulacion-> 1e-3 y 100
% polos=eig(Mat_A)
% tR=log(0.95)/polos(3); %dinámica mas rápida
% tL=log(0.05)/polos(2); %dinámica mas lenta


%CONDICIONES INICIALES
alpha(1) = 1.17;        %El maximo ang inicial
% alpha(1) = 0.1; 
ref = -10;
flag = 0;

Tf=100;h=1e-3;pasos=Tf/h;
p_pp = 0;
tita_pp = 0;
omega(1) = 0; 
p_p(1) = 0; 
u(1) = 0; 
p(1) = 0; 
i = 1;

%METODO LQR
D = [1e2 100 1e4 20];  %Velocidad angular, Posicion angular, Posicion Carrito, Velocidad Carrito 
Q = diag(D);
R = 1000;
Klqr = lqr(Mat_A,Mat_B,Q,R);

while(i<(pasos+1))
    X = [p(i); p_p(i); alpha(i); omega(i)];
        %Ley de control
    u(i) = -Klqr*X+Klqr(1)*ref; color = 'm';              
    
    %Sitema No Lineal

    p_pp = (1/(M+m))*(u(i)-m*long*tita_pp*cos(alpha(i))+m*long*omega(i)^2*sin(alpha(i))-Fricc*p_p(i));
    tita_pp = (1/long)*(g*sin(alpha(i))-p_pp*cos(alpha(i)));
    
    %Integracion por Euler
    p_p(i+1) = p_p(i)+h*p_pp;
    p(i+1) = p(i)+h*p_p(i);
    omega(i+1) = omega(i)+h*tita_pp;
    alpha(i+1) = alpha(i)+h*omega(i);
    
    i = i+1;
    
end

u(i) = u(i-1);
t = 0:pasos; 
t = t*h;

%Graficas 
figure(1);hold on;
subplot(3,2,1); plot(t,omega,color); grid on; title('Velocidad ángulo'); hold on;
subplot(3,2,2); plot(t,alpha,color); grid on; title('Ángulo'); hold on;
subplot(3,2,3); plot(t,p,color); grid on; title('Posición carro'); hold on;
subplot(3,2,4); plot(t,p_p,color); grid on; title('Velocidad carro'); hold on;
subplot(3,1,3); plot(t,u,color); grid on; title('Acción de control'); xlabel('Tiempo en Seg.'); hold on;