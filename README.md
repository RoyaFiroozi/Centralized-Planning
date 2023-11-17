# Centralized-Planning
Optimization-based formation and reconfiguraiton of multi-lane vehicle platoons is performed using Model Predictive Control (MPC).
The collision avoidance constraints among the vehicles (dynamic obstacles) as well as static obstacles on the road are modeled using strong duality theory. This formulation allows for motion planning and obstacle avoidance in tight environments. 
Paper describing the theory can be found [here](https://arxiv.org/abs/2003.08595).

The requirements to run the code are MATLAB, YALMIP, and IPOPT solver for nonlinear optimization.

## Examples
### Obstacle Avoidance Scenario
<img src="https://github.com/RoyaFiroozi/Centralized-Planning/blob/master/Obstacle_Avoidance.gif" width="700" />

### Platoon Reconfiguration Scenario
<img src="https://github.com/RoyaFiroozi/Centralized-Planning/blob/master/Platoon_Reconfiguration.gif" width="700" />

