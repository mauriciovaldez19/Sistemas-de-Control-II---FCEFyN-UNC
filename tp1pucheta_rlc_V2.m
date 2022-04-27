%     Ejercicio 1 Tp Control 2 - Prof: Dr. Pucheta; Alumno: Valdez Benavidez, Mauricio Luciano
%    ------------------------------------------------------------------------------------------------
% 1)  Asignar valores a R=4,7Kohm, L=10uHy, y C=100nF. Obtener simulaciones que
%     permitan estudiar la dinámica del sistema, con una entrada de tensión escalón de 12V,
%     que cada 1ms cambia de signo.
%2)   Asignar valores a R=5,6Kohm, L=10uHy, y C=100nF; repetir lo anterior para comparar el 
%     resultado y verificar la correcta simulación

clc, clear all, close all
% Defino valores de RLC
%R=4.7e3;  L=10e-6; C=100e-9; %para punto 1)
R=5.6e3;  L=10e-6; C=100e-9;  %para punto 2)
% Defino matrices
matA=[[-R/L -1/L];
      [1/C 0]]; 
matB=[1/L ;
        0]; 
matC=[R 0]'; 

% Para obtener las simulaciones, necesito encontrar el valor de los tiempos de muestreo,
% busco los polos de la funcion de transferencia, para buscar la dinamica mas rapida
% como tengo la Matriz A ya definida, obteniendo los autovalores, tengo los polos de la funcion de transferencia

autovalor=eig(matA)

tR=log(0.95)/autovalor(1)
h=tR/3

tL=log(0.05)/autovalor(2)
tsim=tL*3

%simulacion

Vin=-12; Vant=0; ii=1; 
t=[]; paso=tsim/h; %aux=0.001/h;

%Condiciones iniciales
vc(1)=0; il(1)=0; u(1)=12;
tant=1e-3;
X0=[0 0]';x=[0 0]';
Vin=Vin*(-1);
while(ii<=paso+1)
%Variables del sistema lineal
il(ii)=x(1); vc(ii)=x(2);
u(ii)=Vin;
t(ii)=ii*h;
%Sistema lineal
xp=matA*(x-X0)+ matB*u(ii);
x=x+h*xp;
if(t(ii)>tant)
Vin=Vin*(-1); %Con esto varía de 12 a -12 cada 1ms
tant=tant+1e-3;
end
ii=ii+1;
end
subplot(3,1,1);hold on;grid on; %Grafica i, Vc y la entrada
plot(t,il,'r');title('Corriente , t');
subplot(3,1,2);hold on;grid on;
plot(t,vc,'r');title('Vc, t');
subplot(3,1,3);hold on;grid on;
plot(t,u,'r');title('Entrada , t');
