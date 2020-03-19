import numpy as np

class robot:
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
	
	def dynamics(self,current_state):
		dydt=np.matmul(current_state,slef.Jmat)
		return dydt
