clear all; clc;
[x1]=xlsread('Curvas_Medidas_RLC');
t0=x1(:,1);
y=x1(:,3);
t_inic=.029;
[val lugar] =min(abs(t_inic-t0)); y_t1=y(lugar);
t_t1=t0(lugar);
ii=0; % for t_inic=10:15
ii=ii+1;
[val lugar] =min(abs(2*t_inic-t0));
t_2t1=t0(lugar);
y_2t1=y(lugar);
[val lugar] =min(abs(3*t_inic-t0));
t_3t1=t0(lugar);
y_3t1=y(lugar);
K=y(end)/12;
k1=(1/12)*y_t1/K-1;%Afecto el valor del Escalon
k2=(1/12)*y_2t1/K-1;
k3=(1/12)*y_3t1/K-1;
be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));
beta=(k1+alfa2)/(alfa1-alfa2); %(2*k1^3+3*k1*k2+k3-sqrt(be))/(sqrt(be));
T1_ang=-t_t1/log(alfa1);
T2_ang=-t_t1/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;
T1(ii)=T1_ang;
T2(ii)=T2_ang;
T3(ii)=T3_ang;
T3_ang=sum(T3/length(T3));
T2_ang=sum(T2/length(T2));
T1_ang=sum(T1/length(T1));
sys_G_ang=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))
step(12*sys_G_ang,.1)
% step(12*sys_G_ang)
hold on 
plot(t0,y,'r')
% step(sys_G,opt),hold on
% figure