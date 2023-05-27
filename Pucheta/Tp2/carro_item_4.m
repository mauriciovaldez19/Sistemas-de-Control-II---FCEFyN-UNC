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
theta(1) = 0.0; 
ref = -10;
% flag = 0;

%METODO LQR
D = [1e2 100 1e4 20];  %Velocidad angular, Posicion angular, Posicion Carrito, Velocidad Carrito 
Q = diag(D);
R = 1000;
Klqr = lqr(Mat_A,Mat_B,Q,R);

%LQR observador
Qo=1*diag([1 10 1 10]);    Ro=1;
Ao=Mat_A';
Bo=Mat_C';
Co=Mat_B';
Ko=lqr(Ao,Bo, Qo, Ro);

Tf=100;h=1e-3;pasos=Tf/h;
delta_pp = 0; delta_p(1) = 0; delta(1) = 0; 
theta_pp = 0; theta_p(1) = 0; theta(1) = 0;
delta_obs(1)=0;delta_p_obs(1)=0;
u(1) = 0; 

i = 1;
while(i<(pasos+1))
    X = [delta(i); delta_p(i); theta(i); theta_p(i)];
    X_obs = [delta_obs(i); delta_p_obs(i); theta(i); theta_p(i)];
        %Ley de control
    u(i) = -Klqr*X+Klqr(1)*ref;     color = 'r';              
    u_obs(i) = -Ko*X_obs+Ko(1)*ref;      
   
%   SIN OBSERVADOR
    %Sitema No Lineal
    delta_pp = (1/(M+m))*(u(i)-m*long*theta_pp*cos(theta(i))+m*long*theta_p(i)^2*sin(theta(i))-Fricc*delta_p(i));
    theta_pp = (1/long)*(g*sin(theta(i))-delta_pp*cos(theta(i)));
    
    %Integracion por Euler
    delta_p(i+1) = delta_p(i)+h*delta_pp;
    delta(i+1) = delta(i)+h*delta_p(i);
    theta_p(i+1) = theta_p(i)+h*theta_pp;
    theta(i+1) = theta(i)+h*theta_p(i);
    
%   CON OBSERVADOR
    %Sitema No Lineal
    delta_pp_obs = (1/(M+m))*(u_obs(i)-m*long*theta_pp*cos(theta(i))+m*long*theta_p(i)^2*sin(theta(i))-Fricc*delta_p(i));
        
    %Integracion por Euler
    delta_p_obs(i+1) = delta_p(i)+h*delta_pp_obs;
    delta_obs(i+1) = delta(i)+h*delta_p_obs(i);
    

    
    
    i = i+1;
    
end

u(i) = u(i-1);
t = 0:pasos; 
t = t*h;

%Graficas 
figure(1);hold on;
subplot(3,2,1); plot(t,theta_p,color); grid on; title('Velocidad ángulo'); hold on;
subplot(3,2,2); plot(t,theta,color); grid on; title('Ángulo'); hold on;
subplot(3,2,3); plot(t,delta,color); grid on; title('Posición carro'); hold on;
subplot(3,2,3); plot(t,delta_obs,'g'); hold on;
subplot(3,2,4); plot(t,delta_p,color); grid on; title('Velocidad carro'); hold on;
subplot(3,2,4); plot(t,delta_p_obs,'g'); hold on;
subplot(3,1,3); plot(t,u,color); grid on; title('Acción de control'); xlabel('Tiempo en Seg.'); hold on;

