# class for the controllers
# originally written by Prof. Sawyer B Fuller
# converted to python by Daksh Dhingra

import numpy as np
import control.matlab as ctrl

class Controller:
	def __init__(self,dt):
		num = np.array([0.01, 0.01, 0.02])
		den = np.array([0.001, 1, 0])
		sys_tf = ctrl.tf(num, den)
		self.A,self.B,self.C,self.D = ctrl.ssdata(ctrl.c2d(sys_tf,dt,'foh'))
		self.X = np.zeros((self.A.shape[0]))

	def lateral_controller(slef, robot, current_state, desired_state):
		omegabody = current_state[3:6]
		tau_c = -(2*pow(10,-7)) * omegabody
		tau_c = np.maximum(np.minimum(tau_c, robot.max_torque),-robot.max_torque)
		tau_c[2] = 0.0
		return tau_c
	def altitude_controller(self, robot, current_state, desired_state):
		error= desired_state[8]- current_state[8]

		self.X= np.matmul(self.A, self.X)+  self.B.T*error
		self.X =np.asarray(self.X).reshape(-1)
		f_t = np.matmul(self.C, self.X )+ self.D* error
		f_t += 0.8* robot.m* robot.g # feed-forward term
		f_t = np.maximum(np.minimum(f_t, robot.max_f_l),-robot.max_f_l)

		return np.asarray(f_t).reshape(-1)





