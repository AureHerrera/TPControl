function [X]=ModeloRLC(t_etapa, xant, accion)
R=5600; L=10e-6; C=100e-9;
A=[-R/L -1/L;1/C 0]; B=[1/L;0];
h=1e-11;
x= xant;
u= accion;
for i=1:10*t_etapa/h
xp = A*x + B*u;
x = x + h*xp;
end
X=x;