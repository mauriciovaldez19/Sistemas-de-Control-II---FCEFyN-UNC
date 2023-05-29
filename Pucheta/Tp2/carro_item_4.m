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
Mat_C = [1 0 0 0];


%CONDICIONES INICIALES
% alpha(1) = 1.17;        %El maximo ang inicial
alpha(1) = 0.01; 
ref = -10;

h = 1e-3;
tiempo = (100/h);
p_pp = 0;
tita_pp = 0;

omega(1)= 0; p_p(1)= 0; u(1)= 0; p(1)= 0; 

i = 1;

%METODO LQR
D = [1e2 100 1e4 20];  %Velocidad angular, Posicion angular, Posicion Carrito, Velocidad Carrito 
Q = diag(D);
R = 1000;
Klqr = lqr(Mat_A,Mat_B,Q,R);
Glqr=-inv(Mat_C*(eye(4)-inv(Mat_A+Mat_B*Klqr))*Mat_B); %%Es igual a Klqr(1)

%LQR observador
Qo=1*diag([1 10 1 10]);    Ro=1;
Ao=Mat_A';
Bo=Mat_C';
Co=Mat_B';
Ko=lqr(Ao,Bo, Qo, Ro);

x_hat = [0 0 0 0]';

i = 1;

while(i<(tiempo+1))
    
    X = [p(i); p_p(i); alpha(i); omega(i)];
    
    %Ley de control
    u(i) = -Klqr*X+Klqr(1)*ref; color = 'm';              
%   u(i) = -Klqr*x_hat+Klqr(1)*ref; color = 'g'; 
    
    %Sitema No Lineal
    p_pp = (1/(M+m))*(u(i)-m*long*tita_pp*cos(alpha(i))+m*long*omega(i)^2*sin(alpha(i))-Fricc*p_p(i));
    tita_pp = (1/long)*(g*sin(alpha(i))-p_pp*cos(alpha(i)));
    
    %Integracion por Euler
    p_p(i+1) = p_p(i)+h*p_pp;
    p(i+1) = p(i)+h*p_p(i);
    omega(i+1) = omega(i)+h*tita_pp;
    alpha(i+1) = alpha(i)+h*omega(i);
    
%     %Observador
%     y_salO(i) = Mat_C*x_hat;
%     y_sal(i) = Mat_C*X;
%     
%     x_hatp = Mat_A*x_hat+Mat_B*u(i)+Ko*(y_sal(i)-y_salO(i));
%     x_hat = x_hat+h*x_hatp ;
    
    i = i+1;
    
end

u(i) = u(i-1);
t = 0:tiempo; 
t = t*h;

%Graficas 
figure(1);hold on;
subplot(3,2,1); plot(t,omega,color); grid on; title('Velocidad ángulo'); hold on;
subplot(3,2,2); plot(t,alpha,color); grid on; title('Ángulo'); hold on;
subplot(3,2,3); plot(t,p,color); grid on; title('Posición carro'); hold on;
subplot(3,2,4); plot(t,p_p,color); grid on; title('Velocidad carro'); hold on;
subplot(3,1,3); plot(t,u,color); grid on; title('Acción de control'); xlabel('Tiempo en Seg.'); hold on;
