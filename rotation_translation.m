function [A, b] = rotation_translation(x0,theta,h,w)

R = [cos(theta), -sin(theta);
     sin(theta), cos(theta)];
 
b(1) = h/2;
b(2) = w/2;
b(3) = h/2;
b(4) = w/2;
 
A = [R'; -R']; 
b = b' + A*x0;






