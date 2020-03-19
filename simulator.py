import numpy as np
from scipy.integrate import odeint


class Simulator:
	#this class contains have plot, animation, numerical ode solver methods.

	def __init__(self,rob,con):
		self.robot = rob
		self.controller= con

	def test_method(self):
		print(self.robot.J)

	def equatin_solver(self,start_state, desired_state,time, dt):
		current_state= start_state
		time_vector= np.arange(0,time,dt);
		for it,current_step in enumerate(time_vector):
			#tau_c = self.controller.lateral(current_step, current_state, desired_state)
			#f_l = self.controller.altitude(current_step, current_state, desired_state)
			#u= np.vstack((f_l,tau_c))
			current_state=np.reshape(current_state, (3))
			time_interval = np.arange(time_vector[it], time_vector[it]+dt,0.001)
			sol = odeint(self.robot.dynamics, current_state, time_interval)
			current_state = sol[-1]
			
			
	
	
		
