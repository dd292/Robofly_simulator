from simulator import Simulator
from robot import Robot
from controller import Controller
import numpy as np

if __name__ == '__main__':
	start_state = np.array([0, 1, 0])
	desired_state = np.array([0, 0, 0])
	time = 1
	dt = 0.01
	robofly = Robot()
	simple_controller = Controller()
	dsim = Simulator(robofly, simple_controller)
	dsim.equatin_solver(start_state, desired_state, time, dt)
