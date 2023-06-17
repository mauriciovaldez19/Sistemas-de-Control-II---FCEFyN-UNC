% close all; clear all; clc;
%Se declaran las variables

mnom = 1;% masa nominal
m=mnom;color='r';% corro programa con masa nominal
% m=mnom*0.90;color='m'; % correr de nuevo el código de simulación, dibujo y analisis
% m=mnom*1.1;color='k'; % correr de nuevo el código de simulación, dibujo y analisis

b = 0.4;
l = 1;
g = 10;
delta = 180;        %En grados


%linealizacion del sistema
Aa = [         0              1;
      (-g/l)*cosd(delta)   -b/(m/l^2)]

Bb = [    0;
    1/(m*l^2)]

Cc = [1 0]

%autovalores
autoval=eig(Aa)

%comparacion de resultados
[A,B,C,D]=linmod('pendulo_mod_tarea',delta*pi/180)
autoval_2 = eig(A)

%Se verifica la controlabilidad de la matriz A
if (length(A)==rank(ctrb(A,B)))
    
    disp('La matriz A es controlable')

else
    
    disp('La matriz A NO es controlable')
    
end 


%matrices del sistema ampliado
AA=[[A;C] zeros(3,1)]
BB=[B;0]

%autovalores, estabilidad y controlabilidad del sistema ampliado
autoval_AA = eig(AA)
rank(ctrb(AA,BB))

%Se verifica la controlabilidad de la matriz ampliada AA
if (length(AA)==rank(ctrb(AA,BB)))
    
    disp('La matriz AA es controlable')

else
    
    disp('La matriz AA NO es controlable')
    
end 

%diseño de un controlador
%asignacion de polos

p = -4;             %Polo triple
K = acker(AA,BB,[p p p])
k1 = K(1)
k2 = K(2)
k3 = K(3)
eig(AA-BB*K)        %Polos lazo cerrado
tscalc = 7.5/(-p)   %Tiempo de respuesta calculado

%simulacion PID
sim('pendulo_pid_tarea');
figure(1); 
plot(tout,yout,color,'LineWidth',2); grid on; title('Salida');hold on;
legend('m=1','m=0.9','m=1.1');legend('boxoff');


%Plano de fase
figure(2);
plot(yout,velocidad,color,'LineWidth',2); grid on; title('Plano de fases'); hold on;
legend('m=1','m=0.9','m=1.1');legend('boxoff');
%Torque total
figure(3); 
plot(tout,torque,color,'LineWidth',2); grid on;title('Torque');hold on;
legend('m=1','m=0.9','m=1.1');legend('boxoff');
%Accion integral
figure(4);
plot(tout,-accint,color,'LineWidth',2); grid on;title('Accion integral');hold on;
legend('m=1','m=0.9','m=1.1');legend('boxoff');

disp('Maximo valor de salida:')
ymax=max(yout) 

disp('Sobrepaso en %:')
S=(ymax-delta)/delta*100

disp('Error relativo:')
erel=(delta-yout)/delta

disp('Error final, debe ser cero:')
efinal=erel(end) 

disp('Indice elementos con error relativo absoluto menor a 2:')
ind=find(abs(erel)>.02)

disp('Tiempo de establecimiento (ultimo valor del vector):')
tss=tout(ind(end))       

disp('Salida al tiempo ts:')
yte=yout(ind(end))   

disp('Torque final:')
uf=torque(end)

disp('Accion integral final:')
Intf=-accint(end)           


