function [A, b] = rotation_translation(x0, theta, h, w)
% x0 is a 2by1 vector [x,y]' which represents the coordinates of the center of the vehicle
% theta is the heading angle in radian
% h is the total length of the vehicle
% w is the width of the vehicle

R = [cos(theta), -sin(theta); % R is rotation matrix
     sin(theta), cos(theta)];
 
b(1) = h/2;
b(2) = w/2;
b(3) = h/2;
b(4) = w/2;
 
% vehicle road occupied region is modeled as a polytopic set
A = [R'; -R']; % A is polytopic representation of the vehicle
b = b' + A*x0; % b is polytopic representation of the vehicle 






