clear all, close all,clc
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
hold on;
plot(t,u, 'm' );title('Tensión de Entrada, u_t');grid on
