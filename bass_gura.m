%BassGura

clear all
close all
clc

Sn=0.00005; S=0.0154;

a_13=0.4753*Sn*sqrt(2*9.8);
a_32=0.4833*Sn*sqrt(2*9.8);
a_20=0.9142*Sn*sqrt(2*9.8);
Q_10 = 3e-5;
Q_20 = 0.5e-5;

H_20 = ((Q_10 + Q_20)/a_20)^2;
H_30 = H_20 + (Q_10/a_32)^2;
H_10 = H_30 + (Q_10/a_13)^2;

R_13=2*sqrt(abs(H_10-H_30))/a_13;
R_32=2*sqrt(abs(H_30-H_20))/a_32;
R_20=2*sqrt(abs(H_20))/a_20;

A=[-1/(S*R_13) 1/(S*R_13) 0;
    1/(S*R_13) -(1/S)*((1/R_13)+(1/R_32)) 1/(S*R_32);
    0 1/(S*R_32) -(1/S)*((1/R_32)+(1/R_20))];
B=[1/S 0; 0 0;0 1/S];
%C=eye(3);
C = [1 0 0];
D=0;

vp=[-0.05;-0.04;-0.035];
q=[2 ; 0];
Bt=B*q;
Ct=ctrb(A,Bt)
rang_Ct=rank(Ct)
k=acker(A,Bt,vp);
K=q*k;

syscor = ss(A-B*K,B,C,D);
syscortf = zpk(tf(syscor))

Gbg=dcgain(syscortf)

N1bg = 1/Gbg(1)
N2bg = (1-evalfr(syscortf(1),0)*N1bg)/evalfr(syscortf(2),0)
precompbg = [N1bg; N2bg];

new_sys_bf1 = ss(A-B*K, B(:,1), [0 0 0; eye(3)],[1; zeros(3,1)]);
new_sys_bf2 = ss(A-B*K, B(:,2), [0 0 0; eye(3)],[1; zeros(3,1)]);

%[Y,t,x] = step(new_sys_bf1*precompbg, stepDataOptions('StepAmplitude',0.05));

[Y,t,x] = step(syscor*precompbg, stepDataOptions('StepAmplitude',0.05)); % q_1 and q_2 together

figure(1)
plot(t(:),x(:,1))
title('h_1(t)')
