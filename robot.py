# class for the dynamics of robot
# originally written by Prof. Sawyer B Fuller
# converted to python by Daksh Dhingra

import numpy as np

class Robot:
	winglength = .6 * 2.54 / 100 
	l = .013 # length of fly
	h = .0025 # thickness of body
	J = 1.5e-9  # 1d version
	Jmat = np.diag([J, J, .5e-9])  #1.4e-9  # from science paper/solidworks
	b_w = 2.0e-4  # aero drag on wings from wind tunnel tests, Ns/m
	c_w = (h/2 + winglength * 2./3)**2 * b_w  # rot drag coeff around z [Nsm]
	r_w = .007  #p.l  # z-dist from ctr of wings to ctr of mass 
	m = 111e-6  #80e-6  #mass of fly
	g = 9.81  
	max_f_l = 1.5 * m * g  
	max_torque = 2e-6   
	ks = 0 #0.5e-6  # tether stiffness
	force_bias_x = 0  #0.0001  #N
	torque_bias_y = 0  #-0.1e-6  #Nm
	gyro_bias_y = 0  #0.1  # rad/s
	# leave the following three zero for this simulation:
	force_bias_y = 0  #N
	torque_bias_x = 0  #Nm
	gyro_bias_x = 0  # rad/s

#state =[theta, omegabody, posworld, vbody]

	def dynamics(self, current_state, time, u):
		theta = current_state[0:3]
		omegabody = current_state[3:6]
		posworld = current_state[6:9]
		vbody = current_state[9:12]

		f_c = u[0]
		tau_c = u[1:]

		R = self.rot_matrix(theta)
		W = self.omega2thetadot_matrix(theta)
		aerodynamics= self.fly_aero(vbody, omegabody)
		f_d = aerodynamics[0, :].transpose()
		tau_d = aerodynamics[1, :].transpose()

		# Forces in body coords
		f_g = np.matmul(R.transpose(), np.array([0, 0, -self.m*self.g]))
		f_l = np.array([0, 0, f_c])

		f = f_l + f_g + f_d.reshape([3])

		#moment/torques (body coords)
		tau = tau_c + tau_d - self.ks*np.array([[theta[0]], [theta[1]], [0]])
		#fictitious force and torque
		fictitious_f= self.m *np.cross(omegabody,vbody)
		fictitious_tau= np.cross(omegabody, np.matmul(self.Jmat,omegabody))

		#geometric
		xdotworld= np.matmul(R, vbody)
		vdotbody= (1/self.m) * (f-fictitious_f)
		thetadot= np.matmul(W, omegabody)
		vec= (tau.reshape([3])-fictitious_tau)

		omegadotbody= np.linalg.lstsq(self.Jmat, vec)[0]
		qdot = np.array([thetadot, omegadotbody, xdotworld, vdotbody])

		return qdot



		return dydt
	def fly_aero (self, v, omega):
		# calculate stroke-averaged forces due to aerodynamic drag on flapping
		# wings. assumes force is applied at point r_w away from center of mass in
		# z-direction in body-coordinates. v_w is the vel at that point.
		# assumes drag force f_d = -b_w * v_w. this has been tested in a wind
		# tunnel for x- and y-directions but not z-direction
		r_w = np.array([0, 0, self.r_w])
		v_w = v + np.cross(omega, r_w)
		f_d = -self.b_w*v_w
		tau_d= np.cross(r_w, f_d)

		return np.array([[f_d],[tau_d]])






	def rot_matrix(self,theta):
		#convention Rotate about Z by theta3 then rotate about new Y by theta2 and rotate about new X by theta1. Alternatively,
		# rotate about X by theta1 then rotate about inital Y by theta2 and  finally rotate about inital Z by theta3
		cz= np.cos(theta[2])
		cy = np.cos(theta[1])
		cx = np.cos(theta[0])
		sz = np.sin(theta[2])
		sy = np.sin(theta[1])
		sx = np.sin(theta[0])
		R = np.array([[cz*cy,   cz*sy*sx - cx*sz,   sz*sx + cz*cx*sy],
					[cy*sz,   cz*cx + sz*sy*sx,   cx*sz*sy - cz*sx],
					[ -sy,    cy*sx,              cy*cx]])
		return R
	def omega2thetadot_matrix(self,euler_theta):# still trying to figure what it is
		# transform euler angle rates to angular rot rate vector omega
		# thetadot = W * omega

		st1 = np.sin(euler_theta[0])
		ct1 = np.cos(euler_theta[0])
		tt2 = np.tan(euler_theta[1])
		ct2 = np.cos(euler_theta[1])

		#XYZ (airplane) convention
		W = np.array([[1, st1*tt2, ct1*tt2],
			[0, ct1, -st1],
			[0, st1/ct2, ct1/ct2]])
		return W
if __name__ == '__main__':
	start_state = np.array([0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1])
	R1= Robot()
	input= np.array([[1], [1], [1], [1]])
	print(R1.dynamics(start_state, 1, input))