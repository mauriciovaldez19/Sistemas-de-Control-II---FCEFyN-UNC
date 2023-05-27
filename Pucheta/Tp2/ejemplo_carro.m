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

% %tiempo de integracion y tiempo de simulacion-> 2e-3 y 100
% polos=eig(Mat_A)
% tR=log(0.95)/polos(3); %dinámica mas rápida
% tR=tR/10;
% tL=log(0.05)/polos(2); %dinámica mas lenta
% tL=tL*5;

Mat_M = ctrb(Mat_A,Mat_B) %matriz controlabilidad

%CONDICIONES INICIALES
% alpha(1) = 1.17;        %El maximo ang inicial
alpha(1) = 0.0; 
ref = -10;
flag = 0;

h = 1e-3;
tiempo = (100/h);
p_pp = 0;
tita_pp = 0;

omega(1) = 0; 
p_p(1) = 0; 
u(1) = 0; 
p(1) = 0; 

i = 1;


%CONTROLADOR K
auto_val = eig(Mat_A);
pol_caract_A = poly(Mat_A);

Mat_W = [ pol_caract_A(4)  pol_caract_A(3)  pol_caract_A(2)  1;  
          pol_caract_A(3)  pol_caract_A(2)         1         0;
          pol_caract_A(2)         1                0         0;
                1                 0                0         0];
            
Mat_T = Mat_M*Mat_W;
A_controlable = inv(Mat_T)*Mat_A*Mat_T;                 %Verificación de que T esté bien
                                                        %Si la inversa existe es controlable

%METODO LQR
D = [1e2 100 1e4 20]   %Velocidad angular, Posicion angular, Posicion Carrito, Velocidad Carrito 
Q = diag(D);
R = 1000;
[Klqr,Plqr,polos_lqr] = lqr(Mat_A,Mat_B,Q,R);
Glqr = -inv(Mat_C*(eye(4)-inv(Mat_A+Mat_B*Klqr))*Mat_B); %Es igual a Klqr(1)


%OBSERVADOR
Mat_AO = Mat_A';
Mat_BO = Mat_C';
Mat_CO = Mat_B';

Mat_M_Dual = [Mat_BO Mat_AO*Mat_BO Mat_AO^2*Mat_BO Mat_AO^3*Mat_BO]; %Matriz de Controlabilidad del observador
                                                                     

% Ubicacion del Observador
alfa_io = poly(polos_lqr*15);
Mat_TO = Mat_M_Dual*Mat_W;
KO = (fliplr(alfa_io(2:end)-pol_caract_A(2:end))*inv(Mat_TO))';
eig(Mat_AO'-KO*Mat_C);                       %Verifico que todos los polos estén en el semiplano izquierdo

x_hat = [0 0 0 0]';

i = 1;

while(i<(tiempo+1))
    
    X = [p(i); p_p(i); alpha(i); omega(i)];
    
    %Ley de control
%    u(i) = -Klqr*X+Glqr*ref; color = 'm';              
   u(i) = -Klqr*x_hat+Glqr*ref; color = 'r'; 
    
    %Sitema No Lineal
    p_pp = (1/(M+m))*(u(i)-m*long*tita_pp*cos(alpha(i))+m*long*omega(i)^2*sin(alpha(i))-Fricc*p_p(i));
    tita_pp = (1/long)*(g*sin(alpha(i))-p_pp*cos(alpha(i)));
    
    %Integracion por Euler
    p_p(i+1) = p_p(i)+h*p_pp;
    p(i+1) = p(i)+h*p_p(i);
    omega(i+1) = omega(i)+h*tita_pp;
    alpha(i+1) = alpha(i)+h*omega(i);
    
    %Observador
    y_salO(i) = Mat_C*x_hat;
    y_sal(i) = Mat_C*X;
    
    x_hatp = Mat_A*x_hat+Mat_B*u(i)+KO*(y_sal(i)-y_salO(i));
    x_hat = x_hat+h*x_hatp ;
    
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