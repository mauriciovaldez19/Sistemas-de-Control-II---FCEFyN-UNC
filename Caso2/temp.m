clear all;close all; clc;
t_etapa=1e-7;
tF=.6;
u=linspace(0,0,(100e3-1));
ii=0;
for t=0:t_etapa:tF
ii=ii+1;
    u(ii)=0;
    if(t>=0.025)
        u(ii)=12;
    end
%     if(t>=0.1)
%         Tl=7.5e-2;
%     end

     if(t>=0.1501)
         u(ii)=-12;
     end
end
t=0:t_etapa:tF;
plot(t,u)