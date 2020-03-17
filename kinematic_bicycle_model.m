function [x,y,psi,v] = kinematic_bicycle_model(x,y,psi,v,delta_t,l_f,l_r,a,delta_f)

beta = atan(((l_r./(l_f+l_r))*tan(delta_f)));
x = x + delta_t*v*cos(psi+beta); % position in x direction
y = y + delta_t*v*sin(psi+beta); % position in y direction 
% psi = psi + delta_t*(v./l_r)*sin(beta); % heading angle
psi = psi + delta_t*(v*cos(beta)./(l_r+l_f)*tan(delta_f)); % heading angle
v = v + delta_t*a; % vcloselocity (norm of velocity in x and y direction)
end
