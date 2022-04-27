% Adam Elkhanoufi
% CSC 2262
% Lab 6

R1 = 2.56;
R2 = 4.15;
R3 = 3.76;
R4 = 3.25;
guess1 = 70*pi/180;
guess2 = 35*pi/180;
accuracy = 1e-7;
for (t1 = 0 : 30*pi/180 : pi)
 f1 = @(t3,t4) R1*cos(t1) + R4*cos(t4) + R3*cos(t3) - R2;
 f2 = @(t3,t4) R1*sin(t1) + R4*sin(t4) - R3*sin(t3);
 df1dt1 = @(t3,t4) -R3*sin(t3);
 df1dt2 = @(t3,t4) -R4*sin(t4);
 df2dt1 = @(t3,t4) -R3*cos(t3);
 df2dt2 = @(t3,t4)  R4*cos(t4);
 [t3,t4] = Newton2(f1,f2,df1dt1,df1dt2,df2dt1,df2dt2,guess1,guess2,accuracy);
 fprintf('theta1 = %3.0f theta3 = %.5f theta4 = %.5f\n', ...
 t1*180/pi, t3*180/pi, t4*180/pi);
end
% function Newton2
function [t3,t4] = Newton2(f1,f2,df1dt1,df1dt2,df2dt1,df2dt2, ...
 guess1,guess2,accuracy)
t3_new = guess1;
t4_new = guess2;
t3_old = guess1 + 1;
t4_old = guess2 + 1;
while(abs(t3_new - t3_old) >= accuracy || abs(t4_new - t4_old) >= accuracy)
 t3_old = t3_new;
 t4_old = t4_new;
 d = [ f1(t3_old,t4_old)
 f2(t3_old,t4_old) ];
 a = [ df1dt1(t3_old,t4_old) df1dt2(t3_old,t4_old)
 df2dt1(t3_old,t4_old) df2dt2(t3_old,t4_old) ];
 p = a \ d;
 t3_new = t3_old - p(1);
 t4_new = t4_old - p(2);
end
t3 = t3_new;
t4 = t4_new;
end