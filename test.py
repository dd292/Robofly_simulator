from simulator import simulator
from robot import robot
from controller import controller
import numpy as np

if __name__== '__main__':
	start_state= np.array([0,1,0])
	desired_state= np.array([0,0,0])
	time=1
	dt=0.01
	robofly= robot
	simple_controller= controller
	dsim= simulator(robofly,simple_controller)
	dsim.equatin_solver(start_state.T, desired_state.T,time, dt)
