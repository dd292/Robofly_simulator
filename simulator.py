# class for the ode simulator
# originally written by Prof. Sawyer B Fuller
# converted to python by Daksh Dhingra


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import csv
import time

class Simulator:
	#this class contains have plot, animation, numerical ode solver methods.

	def __init__(self,rob,con,write_data):
		self.robot = rob
		self.controller= con
		self.write_data= write_data

	def test_method(self):
		print(self.robot.J)

	def equatin_solver(self,start_state, desired_state,time_var, dt):
		start= time.time()
		current_state= start_state
		time_vector= np.arange(0,time_var,dt);
		states= []
		writer = csv.writer(open("final_states.csv", 'w',newline=''))
		for it,current_step in enumerate(time_vector):
			tau_c = self.controller.lateral_controller(self.robot, current_state, desired_state)
			f_l = self.controller.altitude_controller(self.robot, current_state, desired_state)
			u= np.concatenate((f_l,tau_c))

			time_interval = np.arange(time_vector[it], time_vector[it]+dt,dt/10)
			sol = odeint(self.robot.dynamics, current_state, time_interval, args=(u,))
			current_state = sol[-1]
			states.append(current_state)

			if (self.write_data):
				writer.writerow(current_state)

		end= time.time()
		total_time= end-start
		print(total_time)
		self.plot_results(np.asarray(states),dt)


	def plot_results(self,states,dt):
		final_time = states.shape[0]*dt
		t= np.arange(0, final_time, dt)
		plt.subplot(3,2,1)
		plt.plot(t,states[:,6])
		plt.title('X position vs time')
		plt.ylabel('X Pos')

		plt.subplot(3,2,2)
		plt.plot(t,states[:,7])
		plt.title('Y position vs time')
		plt.ylabel('Y Pos')


		plt.subplot(3,2,3)
		plt.plot(t,states[:,8])
		plt.title('Z position vs time')
		plt.ylabel('Z Pos')


		plt.subplot(3,2,4)
		plt.plot(t,states[:,0])
		plt.title('theta1 vs time')
		plt.ylabel('theta1')

		plt.subplot(3,2,5)
		plt.plot(t,states[:,1])
		plt.title('theta2 vs time')
		plt.ylabel('theta2')

		plt.subplot(3,2,6)
		plt.plot(t,states[:,2])
		plt.title('theta3 vs time')
		plt.ylabel('theta3')

		plt.show()


			
	
	
		
