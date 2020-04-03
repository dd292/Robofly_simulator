from simulator import Simulator
from robot import Robot
from controller import Controller
import numpy as np
import time

if __name__ == '__main__':

	start_state = np.array([.2, -.2, 0,  -1, 0, 1, .04, .04, .01,  .1, -.3, 0])
	desired_state = np.array([0,0,0,  0,0,0, 0,0,.08,  0,0,0])
	time = 10
	dt = 1/10000
	write_data= True
	robofly = Robot()
	simple_controller = Controller(dt)
	dsim = Simulator(robofly, simple_controller,write_data)
	dsim.equatin_solver(start_state, desired_state, time, dt)

