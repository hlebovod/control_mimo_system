function F = myfun(x)

Sn=0.00005; S=0.0154;
a_13=0.4753*Sn*sqrt(2*9.8);
a_32=0.4833*Sn*sqrt(2*9.8);
a_20=0.9142*Sn*sqrt(2*9.8);
Q_10 = 0.00003;
Q_20 = 0.000005;

F = [Q_10^2 - (x(1) - x(2)) * a_13^2;
      a_13 * sqrt(x(1) -x(2)) - a_32 * sqrt(x(2) - x(3));
      Q_20 + a_32 * sqrt(x(2) - x(3)) - a_20 * sqrt(x(3))];

