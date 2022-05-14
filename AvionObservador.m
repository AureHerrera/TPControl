clc;clear all;
w=2;a=0.05;b=5;c=100;
dt=1e-4;tf=20;tiempo=(tf/dt); t=0:dt:tf;
alfa=0:dt:tf; fi=0:dt:tf; fi_p=0:dt:tf; h=0:dt:tf;
y_sal=0:dt:tf; u=linspace(0,dt,tiempo+1);
%condiciones iniciales
alfa(1)=0.2; fi(1)=0.3; h(1)=500; i=1; color='r';
ref=100; indice=0;
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[alfa(i); fi(i); fi_p(i); h(i)]
Mat_A=[-a a 0 0;0 0 1 0; w^2 -w^2 0 0; c 0 0 0];
Mat_B=[0; 0; b*w^2; 0];
Mat_C=[0 0 0 1]; %La salida es altura
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B];%Matriz Controlabilidad
%Cálculo del controlador por asignación de polos
auto_val=eig(Mat_A);
c_ai=conv(conv(conv([1 -auto_val(1)],[1 -auto_val(2)]),[1 -auto_val(3)]),[1 -auto_val(4)]);
Mat_W=[c_ai(4) c_ai(3) c_ai(2) 1;c_ai(3) c_ai(2) 1 0;c_ai(2) 1 0 0;1 0 0 0];
% Mat_W=[c_ai(3) c_ai(2) c_ai(1) 1;c_ai(2) c_ai(1) 1 0;c_ai(1) 1 0 0;1 0 0 0];
Mat_T=Mat_M*Mat_W;
B_controlable=inv(Mat_T)*Mat_B;
A_controlable=inv(Mat_T)*Mat_A*Mat_T; %Verificación de que T esté bien
%CONTROLADOR Ubicación de los polos de lazo cerrado en mui :
% mui(1)=-50;mui(2)=-30; mui(3)=-1 + i;mui(4)=conj(mui(3));
mui(1)=-50;mui(2)=-30; mui(3)=-1+i;mui(4)=conj(mui(3));
alfa_i=conv(conv(conv([1 -mui(3)],[1 -mui(4)]),[1 -mui(2)]),[1 -mui(1)]);
K=fliplr(alfa_i(2:5)-c_ai(2:5))*inv(Mat_T);
eig(Mat_A-Mat_B*K)
Mat_A_O=Mat_A';
Mat_B_O=Mat_C';
Mat_M_Dual=[Mat_B_O Mat_A_O*Mat_B_O Mat_A_O^2*Mat_B_O Mat_A_O^3*Mat_B_O];%Matriz Controlabilidad
alfaO_i=alfa_i;
% Ubicacion del Observador
% Algunas veces más rápido que el controlador
mui_o=real(mui)*20;
alfaO_i=conv(conv(conv([1 -mui_o(3)],[1 -mui_o(4)]),[1 -mui_o(2)]),[1 -mui_o(1)]);
Mat_T_O=Mat_M_Dual*Mat_W;
Ko=(fliplr(alfaO_i(2:end)-c_ai(2:end))*inv(Mat_T_O))';
eig(Mat_A_O'-Ko*Mat_C) %Verifico que todos los polos estén en el semiplano izquierdo
x_hat=[0;0;0;0]; %Inicializo el Observador
while(i<(tiempo+1))
estado=[alfa(i); fi(i); fi_p(i); h(i)];
% u(i)=-K*estado; color='b'; %Sin Observador
u(i)=-K*x_hat; color='r'; %Con Observador
% Ecuaciones diferenciales
alfa_p = a*(fi(i)-alfa(i));
fi_pp = -(w^2)*(fi(i)-alfa(i)-b*u(i));
h_p = c*alfa(i);
% Integraciones por Euler
alfa(i+1) = alfa(i)+alfa_p*dt;
h(i+1) = h(i)+h_p*dt; 
fi_p(i+1) = fi_p(i)+fi_pp*dt;
fi(i+1) = fi(i)+fi_p(i)*dt;
% y_sal(i)=Mat_C*estado;
%________OBSERVADOR__________
y_sal_O(i)=Mat_C*x_hat;
y_sal(i)=Mat_C*estado;
x_hatp=Mat_A*x_hat+Mat_B*u(i)+Ko*(y_sal(i)-y_sal_O(i));
x_hat=x_hat+dt*x_hatp;
i=i+1;
end
figure(1);hold on; %t=1:i;t=t*dt;
subplot(3,2,1);plot(t,h,color);grid on;title('altura(h)');hold on;
legend('Sin Observador','Con Observador');legend('boxoff');
subplot(3,2,2);plot(t,fi_p,color);grid on;title('Velocidad Ángulo');hold on;
subplot(3,2,3);plot(t,fi,color);grid on;title('direccion de vuelo(fi)');hold on;
subplot(3,2,4);plot(t,alfa,color);grid on; title('direccion de vuelo(alfa)');hold on;
subplot(3,1,3);plot(t,u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;
% figure(2);hold on;
% subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
% subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;
% legend('Sin Observador','Con Observador');legend('boxoff');