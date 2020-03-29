function [x,y,psi,v] = kinematic_bicycle_model(x,y,psi,v,delta_t,l_f,l_r,a,delta_f)
% x is position in longitudinal direction
% y is position in lateral direction
% psi is heading angle
% v is velocity (norm of velocity in x and y directions)
% delta_t is sampling time
% l_f is the length of the car from center of gravity to the front end
% l_r is the length of the car from center of gravity to the rear end
% a is acceleration which is control input
% delta_f is steering angle which is control input 

beta = atan(((l_r./(l_f+l_r))*tan(delta_f)));
x = x + delta_t*v*cos(psi+beta); 
y = y + delta_t*v*sin(psi+beta); 
% psi = psi + delta_t*(v./l_r)*sin(beta); 
psi = psi + delta_t*(v*cos(beta)./(l_r+l_f)*tan(delta_f)); 
v = v + delta_t*a; 
end
