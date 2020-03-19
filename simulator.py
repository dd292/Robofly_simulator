import numpy as np
from scipy.integrate import odeint


class simulator:
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
			sol = odeint(self.robot.dynamics, current_state, np.linspace[time_vector[it],0.001,time_vector[it]+dt])
			print(sol)
			current_state=sol
			
			
	
	
		
