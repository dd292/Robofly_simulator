import matplotlib.pyplot as plt
import numpy as np
import csv

def plot_results(states,dt):
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
if __name__ == '__main__':
    dt= 1/10000
    with open('final_states.csv',  newline='') as f:
        data = list(csv.reader(f))
        all_data=[]
        for row in data:
            s=[]
            for i in row:
                s.append(float(i))
            all_data.append(s)

    all_data= np.asarray(all_data)
    plot_results(all_data,dt)



