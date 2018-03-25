
clc;
clear all;

%Right Arm

l1 = 57;
l2 = 45;
l3 = 75;

%initial angle
px_ini = 0;
py_ini = -177;
pz_ini = -12;

theta1_ini = 90;
theta3_ini = 0;
theta5_ini = 0;


%offset
offset_1 = -90;
offset_3 = 90;
offset_5 = 90;


px = 109;
py = -2;
pz = -25;
    
if (px < 0 && pz > 0) || (px > 0 && pz < 0)
    theta1_1 = atand(pz/px);
    theta1_2 = atand(pz/px) + 180;
elseif (px < 0 && pz < 0) || (px > 0 && pz > 0)
    theta1_1 = atand(pz/px);
    theta1_2 = atand(pz/px) - 180;
end
    

if ((-90 <= theta1_1) || (theta1_1 <= 150)) && ((-90 > theta1_2) || (theta1_2 > 150))
    theta1 = theta1_1;
elseif ((-90 > theta1_1) || (theta1_1 > 150)) && ((-90 <= theta1_2) || (theta1_2 <= 150))
    theta1 = theta1_2;
elseif ((-90 <= theta1_1) || (theta1_1 <= 150)) && ((-90 <= theta1_2) || (theta1_2 <= 150))
    if (theta1_1 > theta1_2) && (theta1_2 > 0)
        theta1 = theta1_2;
    elseif (theta1_1 > theta1_2) && (theta1_2 < 0)
        theta1 = theta1_1;
    elseif (theta1_1 < theta1_2) && (theta1_1 > 0)
        theta1 = theta1_1;
    elseif (theta1_1 < theta1_2) && (theta1_1 < 0)
        theta1 = theta1_2;
    end
end

K = (((px/(3*cosd(theta1)))-4)^2+(py/3+19)^2-850)/750;

theta5_1 = atan2d(K, sqrt(1-K^2));
theta5_2 = atan2d(K, -sqrt(1-K^2));

   
    
    
    if ((-120 <= theta5_1) || (theta5_1 <= 100)) && ((-120 > theta5_2) || (theta5_2 > 100))
        theta5 = theta5_1; %theta5_1 in
    elseif ((-120 > theta5_1) || (theta5_1 > 100)) && ((-120 <= theta5_2) || (theta5_2 <= 100))
        theta5 = theta5_2;
    elseif ((-120 <= theta5_1) || (theta5_1 <= 100)) && ((-120 <= theta5_2) || (theta5_2 <= 100))
        if (theta5_1 > theta5_2) && (theta5_2 > 0)
            theta5 = theta5_2;
        elseif (theta5_1 > theta5_2) && (theta5_2 < 0)
            theta5 = theta5_1;
        elseif (theta5_1 < theta5_2) && (theta5_1 > 0)
            theta5 = theta5_1;
        elseif (theta5_1 < theta5_2) && (theta5_1 < 0)
            theta5 = theta5_2; 
        else
            theta5 = theta5_1;
        end
    end
    
    
   

    
%     a = 25*sind(theta5) + 15;
%     c = 25*sind(theta5);
%     d = px/3*cosd(theta1) - 4;
%     e = 25*cosd(theta5);
%     f = -25*sind(theta5);
%     g = py + 19;
%     
%     theta3 = atan2d(a*g - d*e, d*f - c*g);

    %case 4
    
%     a = (px / (3*cosd(theta1))) - 4;
%     
%     theta3_1 = atan2d(sqrt( (25 * sind(theta5) + 15)^2 + (25 * cosd(theta5))^2 - a^2), a) + atan2d( (25 * cosd(theta5)), (25 * sind(theta5) + 15));
%     theta3_2 = atan2d(-sqrt( (25 * sind(theta5) + 15)^2 + (25 * cosd(theta5))^2 - a^2), a) + atan2d( (25 * cosd(theta5)), (25 * sind(theta5) + 15));
    
    b= py + 57;
    
    theta3_1 = atan2d(sqrt( (75 * cosd(theta5))^2 + (-sind(theta5) - 45)^2 - b^2), b) + atan2d( (-sind(theta5) - 45), (75 * cosd(theta5)));
    theta3_2 = atan2d(-sqrt( (75 * cosd(theta5))^2 + (-sind(theta5) - 45)^2 - b^2), b) + atan2d( (-sind(theta5) - 45), (75 * cosd(theta5)));

%     c = (pz / (3*sind(theta1))) - 4;
%      
%     theta3_1 = atan2d(sqrt( (25 * sind(theta5) + 15)^2 + (25 * cosd(theta5))^2 - c^2), c) + atan2d( (25 * cosd(theta5)), (25 * sind(theta5) + 15));
%     theta3_2 = atan2d(-sqrt( (25 * sind(theta5) + 15)^2 + (25 * cosd(theta5))^2 - c^2), c) + atan2d( (25 * cosd(theta5)), (25 * sind(theta5) + 15)); 

    if ((-100 <= theta3_1) || (theta3_1 <= 45)) && ((-100 > theta3_2) || (theta3_2 > 45))
        theta3 = theta3_1;
    elseif ((-100 > theta3_1) || (theta3_1 > 45)) && ((-100 <= theta3_2) || (theta3_2 <= 45))
        theta3 = theta3_2;
    elseif ((-100 <= theta3_1) || (theta3_1 <= 45)) && ((-100 <= theta3_2) || (theta3_2<= 45))
        if (theta3_1 > theta3_2) && (theta3_2 > 0)
            theta3 = theta3_2;
        elseif (theta3_1 > theta3_2) && (theta3_2 < 0)
            theta3 = theta3_1;
        elseif (theta3_1 < theta3_2) && (theta3_1 > 0)
            theta3 = theta3_1;
        elseif (theta3_1 < theta3_2) && (theta3_1 < 0)
            theta3 = theta3_2; 
        else
            theta3 = theta3_1;
        end
    end
    
    %for motor angle offset
    theta1_final = theta1 - offset_1
    theta3_final = theta3 - offset_3
    theta5_final = theta5 - offset_5
    
% px_ini = 0;
% py_ini = -177;
% pz_ini = -12;
%rmotion
% theta1 = 77
%theta3 = -96
%theta5 = -28

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%Left Arm

%initial angle
px_ini_l = 0;
py_ini_l = 177;
pz_ini_l = -12;

theta2_ini_l = 0;
theta4_ini_l = 0;
theta6_ini_l = 0;


%offset
offset_2 = 90;
offset_4 = -90;
offset_6 = -90;




px_l = 109;
py_l = 2;
pz_l = -25;
    
if (px_l < 0 && pz_l > 0) || (px_l > 0 && pz_l < 0)
    theta2_1 = atand(-pz_l/px_l);
    theta2_2 = atand(-pz_l/px_l) - 180;
elseif (px_l < 0 && pz_l < 0) || (px_l > 0 && pz_l > 0)
    theta2_1 = atand(-pz_l/px_l);
    theta2_2 = atand(-pz_l/px_l) + 180;
end


    if ((-150 <= theta2_1) || (theta2_1 <= 90)) && ((-150 > theta2_2) || (theta2_2 > 90))
        theta2 = theta2_1;
    elseif ((-150 > theta2_1) || (theta2_1 > 90)) && ((-150 <= theta2_2) || (theta2_2 <= 90))
        theta2 = theta2_2;
    elseif ((-150 <= theta2_1) || (theta2_1 <= 90)) && ((-150 <= theta2_2) || (theta2_2 <= 90))
        if (theta2_1 > theta2_2) && (theta2_2 > 0)
            theta2 = theta2_2;
        elseif (theta2_1 > theta2_2) && (theta2_2 < 0)
            theta2 = theta2_1;
        elseif (theta2_1 < theta2_2) && (theta2_1 > 0)
            theta2 = theta2_1;
        elseif (theta2_1 < theta2_2) && (theta2_1 < 0)
            theta2 = theta2_2;
        end
    end
    
    
    
    
M = -(((px_l/(3*cosd(theta2)))-4)^2+(py_l/3-19)^2-850)/750;

theta6_1 = atan2d(M, sqrt(1-M^2));
theta6_2 = atan2d(M, -sqrt(1-M^2));

   
    
    
    if ((-100 <= theta6_1) || (theta6_1 <= 120)) && ((-100 > theta6_2) || (theta6_2 > 120))
        theta6 = theta6_1; %theta6_1 in
    elseif ((-100 > theta6_1) || (theta6_1 > 120)) && ((-100 <= theta6_2) || (theta6_2 <= 120))
        theta6 = theta6_2;
    elseif ((-100 <= theta6_1) || (theta6_1 <= 120)) && ((-100 <= theta6_2) || (theta6_2 <= 120))
        if (theta6_1 > theta6_2) && (theta6_2 > 0)
            theta6 = theta6_2;
        elseif (theta6_1 > theta6_2) && (theta6_2 < 0)
            theta6 = theta6_1;
        elseif (theta6_1 < theta6_2) && (theta6_1 > 0)
            theta6 = theta6_1;
        elseif (theta6_1 < theta6_2) && (theta6_1 < 0)
            theta6 = theta6_2; 
        else
            theta6 = theta6_1;
        end
    end
    
    
   

    b_l= py_l - 57;
    
    theta4_1 = atan2d(sqrt( (-75 * cosd(theta6))^2 + (75*sind(theta6) - 45)^2 - b_l^2), b_l) + atan2d( (75*sind(theta6) - 45), (-75 * cosd(theta6)));
    theta4_2 = atan2d(-sqrt( (-75 * cosd(theta6))^2 + (75*sind(theta6) - 45)^2 - b_l^2), b_l) + atan2d( (75*sind(theta6) - 45), (-75 * cosd(theta6)));

    if ((-45 <= theta4_1) || (theta4_1 <= 100)) && ((-45 > theta4_2) || (theta4_2 > 100))
        theta4 = theta4_1;
    elseif ((-45 > theta4_1) || (theta4_1 > 100)) && ((-45 <= theta4_2) || (theta4_2 <= 100))
        theta4 = theta4_2;
    elseif ((-45 <= theta4_1) || (theta4_1 <= 100)) && ((-45 <= theta4_2) || (theta4_2<= 100))
        if (theta4_1 > theta4_2) && (theta4_2 > 0)
            theta4 = theta4_2;
        elseif (theta4_1 > theta4_2) && (theta4_2 < 0)
            theta4 = theta4_1;
        elseif (theta4_1 < theta4_2) && (theta4_1 > 0)
            theta4 = theta4_1;
        elseif (theta4_1 < theta4_2) && (theta4_1 < 0)
            theta4 = theta4_2; 
        else
            theta4 = theta4_1;
        end
    end
    
    %for motor angle offset
    theta2_final = theta2 - offset_2
    theta4_final = theta4 - offset_4
    theta6_final = theta6 - offset_6
    
% px_ini = 0;
% py_ini = -177;
% pz_ini = -12;
%rmotion
% theta1 = 77
%theta3 = -96
%theta5 = -28

%draw ball

r = 20;
n = 50;
phi = (-n:1:n)'/n*pi;
theta = (-n:1:n)/n*2*pi;
u = r*cos(phi);
x = sqrt(r^2-u.^2)*cos(theta);
y = sqrt(r^2-u.^2)*sin(theta);
z = repmat(u,1,length(phi));
figure(); 

% xlim([-200 200]); 
% ylim([-200 200]);
% zlim([-200 200]); 

surf(x + 109,y,z-25);
view(45, 10);

hold on;

% x_1 = sqrt((r/2)^2-u.^2)*cos(theta);
% y_1 = sqrt((r/2)^2-u.^2)*sin(theta);
% z_1 = repmat(u,1,length(phi));
% surf(x,y,z+5);
% hold on
% grid;

leg_x = [0, 0, 0, 0];
leg_y = [0, 0, -118, 118];
leg_z = [0, -100, -100, -100];


plot3(leg_x, leg_y, leg_z, '-');


title('project 2')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')




aviobj = VideoWriter('test.avi');
open(aviobj);

    
for i = 0: 0.5: 90
    

    
    
    %for forward kinematics
    theta1 = ((i/90) * theta1_final) + offset_1;
    theta2 = ((i/90) * theta2_final) + offset_2;

    theta3 = ((i/90) * theta3_final) + offset_3;
    theta4 = ((i/90) * theta4_final) + offset_4;    

    theta5 = ((i/90) * theta5_final) + offset_5;
    theta6 = ((i/90) * theta6_final) + offset_6;
    
    x3 = 12 * cosd(theta1);
    y3 = -57;
    z3 = 12 * sind(theta1);
    
    x4 = 12 * cosd(theta2);
    y4 = 57;
    z4 = -12 * sind(theta2);
    
    x5 = 12 * cosd(theta1) + 45*cosd(theta1)*cosd(theta3);    
    y5 = -45 * sind(theta3) - 57;
    z5 = 12 * sind(theta1) + 45*cosd(theta3)*sind(theta1);
    
    x6 = 12 * cosd(theta2) + 45*cosd(theta2)*cosd(theta4);    
    y6 = 57 - 45 * sind(theta4);
    z6 = -12 * sind(theta2) - 45*cosd(theta4)*sind(theta2);
    
    x = 3*cosd(theta1)*(25*sind(theta3 + theta5) + 15*cosd(theta3) + 4);
    y = 75*cosd(theta3 + theta5) - 45*sind(theta3) -57;
    z = 3*sind(theta1)*(25*sind(theta3 + theta5) + 15*cosd(theta3) + 4);
    
    x_l = 3*cosd(theta2)*(-25*sind(theta4 + theta6) + 15*cosd(theta4) + 4);
    y_l = -75*cosd(theta4 + theta6) - 45*sind(theta4) + 57;
    z_l = 3*sind(theta2)*(25*sind(theta4 + theta6) - 15*cosd(theta4) - 4);
    

    rx = [0, 0, x3, x5, x];
    ry = [0, -57, y3, y5, y];
    rz = [0, 0, z3, z5, z];
  
    rx_l = [0, 0, x4, x6, x_l];
    ry_l = [0, 57, y4, y6, y_l];
    rz_l = [0, 0, z4, z6, z_l];
    
    axis([-200, 200, -200, 200, -200, 200], 'vis3d'); 
    hold on;
    rightarm = plot3(rx, ry, rz, '-');
    righthand = plot3(x, y, z, '*');
    leftarm = plot3(rx_l, ry_l, rz_l, '-');
    lefthand = plot3(x_l, y_l, z_l, '*');
    
    mo = getframe;
    writeVideo(aviobj, mo);
    
    delete(rightarm);
    delete(leftarm);
    delete(righthand);
    delete(lefthand);

end
    
plot3(rx, ry, rz, '-');
plot3(x, y, z, '*');
plot3(rx_l, ry_l, rz_l, '-');
plot3(x_l, y_l, z_l, '*'); 
axis([-200, 200, -200, 200, -200, 200], 'vis3d'); 
close(aviobj);
    
    
