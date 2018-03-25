%hw3_263a
clear;
clc;

syms c1 c2 c3 c4 s1 s2 s3 s4 c12 s12 c123 s123;




alpha0 = 0;
alpha1 = 0;
alpha2 = 45;
alpha3 = 0;

a0 = 0;
a1 = 1;
a2 = 0;
a3 = sqrt(2);

d1 = 0;
d2 = 0;
d3 = sqrt(2);
d4 = 0;


%  c1 = cosd(0);
%  s1 = sind(0);
%  
%  c2 = cosd(90);
%  s2 = sind(90);
%  
%  c3 = cosd(-90);
%  s3 = sind(-90);
%  
%  c4 = cosd(0);
%  s4 = sind(0);

 
T_01 = [c1, -s1, 0, a0; s1*cosd(alpha0), c1*cosd(alpha0), -sind(alpha0), -sind(alpha0)*d1; s1*sind(alpha0), c1*sind(alpha0), cosd(alpha0), cosd(alpha0)*d1; 0 0 0 1];
 
T_12 = [c2, -s2, 0, a1; s2*cosd(alpha1), c2*cosd(alpha1), -sind(alpha1), -sind(alpha1)*d2; s2*sind(alpha1), c2*sind(alpha1), cosd(alpha1), cosd(alpha1)*d2; 0 0 0 1];

T_23 = [c3, -s3, 0, a2; s3*cosd(alpha2), c3*cosd(alpha2), -sind(alpha2), -sind(alpha2)*d3; s3*sind(alpha2), c3*sind(alpha2), cosd(alpha2), cosd(alpha2)*d3; 0 0 0 1];

T_34 = [c4, -s4, 0, a3; s4*cosd(alpha3), c4*cosd(alpha3), -sind(alpha3), -sind(alpha3)*d4; s4*sind(alpha3), c4*sind(alpha3), cosd(alpha3), cosd(alpha3)*d4; 0 0 0 1];


T_02 = T_01 * T_12;

T_02_1 = [c12, -s12, 0, c1; s12, c12, 0, s1; 0, 0, 1, 0; 0, 0, 0, 1];

T_23_1 = [c3, -s3, 0, 0; 0.7071*s3, 0.7071*c3, -0.7071, -1; 0.7071*s3, 0.7071*c3, 0.7071, 1; 0, 0, 0, 1];

T_03 = T_02 * T_23;
T_03_1 = T_02_1 * T_23_1;

T_34_1 = [ c4, -s4, 0, sqrt(2); s4,  c4, 0, 0; 0, 0, 1,0; 0,0, 0,1];

T_04_1 = T_03_1 * T_34_1;
T_04 = T_01 * T_12 * T_23 * T_34;
T_04 = T_02_1 * T_23_1* T_34_1;


atan2d(-sqrt(1+4-1.4641), 1.21) + atan2d(2, 1);

atan2d(sqrt((092^2)+ (1.115)^2 - (1.1^2)), 1.1) + atan2d(-1.115, 0.92);

atan2d(sqrt((2.048^2)+ (0.389)^2 - (1.1^2)), 1.1) + atan2d(0.3891, 2.048);

atan2d(sqrt((1.6^2)+ (0.94)^2 - (1.1^2)), 1.1) + atan2d(-0.94, 1.6);

px = 1.1;
py = 1.5;
pz = 1.707;

theta3_1 = atan2d(pz - 1, sqrt(1-(pz-1)^2));
theta3_2 = atan2d(pz - 1, -sqrt(1-(pz-1)^2));

a_1 = px^2 + py^2 -3 + sind(theta3_1) - (cosd(theta3_1)^2);
a_2 = px^2 + py^2 -3 + sind(theta3_2) - (cosd(theta3_2)^2);

theta2_1_1 = atan2d(sqrt((2 * sqrt(2) *cosd(theta3_1) )^2 + (- 2 * sind(theta3_1) + 2 )^2 - a_1^2), a_1) + atan2d((-2 * sind(theta3_1) + 2 ), (2*sqrt(2)*cosd(theta3_1)))
theta2_1_2 = atan2d(-sqrt((2*sqrt(2)*cosd(theta3_1))^2 + (- 2 * sind(theta3_1) + 2 )^2 - a_1^2), a_1) + atan2d((-2 * sind(theta3_1) + 2 ), (2*sqrt(2)*cosd(theta3_1)))

theta2_2_1 = atan2d(sqrt((2*sqrt(2)*cosd(theta3_2))^2 + (- 2 * sind(theta3_2) + 2 )^2 - a_1^2), a_1) + atan2d((-2 * sind(theta3_2) + 2 ), (2*sqrt(2)*cosd(theta3_2)))
theta2_2_2 = atan2d(-sqrt((2*sqrt(2)*cosd(theta3_2))^2 + (- 2 * sind(theta3_2) + 2 )^2 - a_1^2), a_1) + atan2d((-2 * sind(theta3_2) + 2 ), (2*sqrt(2)*cosd(theta3_2)))

b_1 = 1+ sind(theta2_1_1) * (1 - cosd(theta3_1)) + sqrt(2) * cosd(theta3_1) * cosd(theta2_1_1);
b_2 = 1+ sind(theta2_1_2) * (1 - cosd(theta3_1)) + sqrt(2) * cosd(theta3_1) * cosd(theta2_1_2);
b_3 = 1+ sind(theta2_2_1) * (1 - cosd(theta3_2)) + sqrt(2) * cosd(theta3_2) * cosd(theta2_2_1);
b_4 = 1+ sind(theta2_2_2) * (1 - cosd(theta3_2)) + sqrt(2) * cosd(theta3_2) * cosd(theta2_2_2);

c_1 = (1 - sind(theta3_1)) * cosd(theta2_1_1) + sqrt(2) * cosd(theta3_1) * sind(theta2_1_1);
c_2 = (1 - sind(theta3_1)) * cosd(theta2_1_2) + sqrt(2) * cosd(theta3_1) * sind(theta2_1_2);
c_3 = (1 - sind(theta3_2)) * cosd(theta2_2_1) + sqrt(2) * cosd(theta3_2) * sind(theta2_2_1);
c_4 = (1 - sind(theta3_2)) * cosd(theta2_2_2) + sqrt(2) * cosd(theta3_2) * sind(theta2_2_2);

d = px;

theta1_1_1_1 = atan2d(sqrt(b_1^2 + c_1^2 - d^2), d) + atan2d(c_1, b_1)
theta1_1_1_2 = atan2d(-sqrt(b_1^2 + c_1^2 - d^2), d) + atan2d(c_1, b_1)

theta1_1_2_1 = atan2d(sqrt(b_2^2 + c_2^2 - d^2), d) + atan2d(c_2, b_2)
theta1_1_2_2 = atan2d(-sqrt(b_2^2 + c_2^2 - d^2), d) + atan2d(c_2, b_2)

%theta1_1_3_1 = atan2d(sqrt(b_3^2 + c_3^2 - d^2), d) + atan2d(c_3, b_3)
%theta1_1_3_2 = atan2d(-sqrt(b_3^2 + c_3^2 - d^2), d) + atan2d(c_3, b_3)

theta1_1_4_1 = atan2d(sqrt(b_4^2 + c_4^2 - d^2), d) + atan2d(c_4, b_4)
theta1_1_4_2 = atan2d(-sqrt(b_4^2 + c_4^2 - d^2), d) + atan2d(c_4, b_4)





