clear all
close all
clc

Sn=0.00005; S=0.0154;
%H_10=0.27474; H_20=0.0299; H_30=0.1368; H_00=0;
%H_10 = 0.19; H_20 = 0.107; H_30 =  0.0299;
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
C=eye(3);
sys=ss(A,B,C,zeros(3,2));
sys_org = ss(A, B, [1 0 0], 0);
sys1=ss(sys.a, sys.b(:,1), [0 0 0; sys.c], [1; sys.d(:,1)]); % q_1
sys2=ss(sys.a, sys.b(:,2), [0 0 0; sys.c], [1; sys.d(:,1)]); % q_2
%sys3=ss(sys.a, sys.b, [0 0 0; sys.c], [1 1; sys.d]); % q_1 with q_2

[Y,t,x]=step(sys1, stepDataOptions('StepAmplitude',5e-5));

% u=ones(3001,2);
% u(:,1)=u(:,1)*5e-5;
% u(:,2)=u(:,2)*2e-5;
% t=0:1:3000;
% [Y,t,x] = lsim(sys,u,t);

subplot(2,2,1)
plot(t(:),Y(:,1))
%plot(t(:), u(:,:))
title('q_2(t)')

subplot(2,2,2)
plot(t(:),x(:,1))
title('h_1(t)')

subplot(2,2,3)
plot(t(:),x(:,2))
title('h_3(t)')

subplot(2,2,4)
plot(t(:),x(:,3))
title('h_2(t)')

%Com=ctrb(sysobs);

tf_norm=tf(sys_org);
ft = (zpk(tf_norm))
