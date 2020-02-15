

import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random



class ParticleBox:
    # Bounds - size of the box [xmin,xmax,ymin,ymax]
    #init_state: [x1,y1,vx1,vy1]
    def __init__(self, init_state = [0,0,0,0],
               bounds = [-4,4,-4,4],
               size = 0.05):
        self.init_state = np.asarray(init_state, dtype = float)
        self.size = size
        self.state = self.init_state.copy()
        self.time_elapsed = 0
        self.bounds = bounds
        self.angle = 0
        
        
    
    def step(self, dt):
        #step once by dt seconds
        self.time_elapsed += dt
        #update positions
        
        self.state[:, :2] += dt*self.state[:, 2:]
        
        #check for crossing boundary
        
        crossed_x1 = (self.state[:, 0] < (self.bounds[0] + self.size))
        crossed_x2 = (self.state[:, 0] > (self.bounds[1] - self.size))
        crossed_y1 = (self.state[:, 1] < (self.bounds[2] + self.size))
        crossed_y2 = (self.state[:, 1] > (self.bounds[3] - self.size))
        
        self.state[crossed_x1, 0] = self.bounds[0] + self.size
        self.state[crossed_x2, 0] = self.bounds[1] - self.size

        self.state[crossed_y1, 1] = self.bounds[2] + self.size
        self.state[crossed_y2, 1] = self.bounds[3] - self.size

        self.state[crossed_x1 | crossed_x2, 2] *= -1
        self.state[crossed_y1 | crossed_y2, 3] *= -1
       
        
        #add force
        
        self.state[:, 3] += -5*dt
      
    def cstep(self, dt):
        self.time_elapsed += dt
        self.angle += 0.01*dt
        self.state[:, 0] = 3*np.cos(self.angle)
        
        self.state[:, 1] = 3*np.sin(self.angle)
        print (self.angle)
        print (self.state)
        return self.state
        
    
    def cforce(self, dt):
        self.time_elapsed += dt

        if (self.time_elapsed >= 6):
            return
        
        self.state[:, 2] += -(1/np.sqrt(((self.state[:, 0])**2)+((self.state[:, 1])**2)))*dt*self.state[:, 0]
        self.state[:, 3] += -(1/np.sqrt(((self.state[:, 0])**2)+((self.state[:, 1])**2)))*dt*self.state[:, 1]
        
        self.state[:, 0] += dt*self.state[:, 2]
        self.state[:, 1] += dt*self.state[:, 3]
        #print(self.state)
        



    def chase(self, dt):
        self.time_elapsed += dt
        #caught:
        caughtx = np.abs(self.state[0,0]-self.state[1,0]) <= 2*self.size
        caughty = np.abs(self.state[0,1]-self.state[1,1]) <= 2*self.size
        if caughtx and caughty:
            return

        #circular motion:
        """
        self.state[0, 2] += -(1/np.sqrt(((self.state[0, 0])**2)+((self.state[0, 1])**2)))*dt*self.state[0, 0]
        self.state[0, 3] += -(1/np.sqrt(((self.state[0, 0])**2)+((self.state[0, 1])**2)))*dt*self.state[0, 1]
        self.state[0, 0] += dt*self.state[0, 2]
        self.state[0, 1] += dt*self.state[0, 3]
        """
        
        crossed_x1 = (self.state[:, 0] < (self.bounds[0] + self.size))
        crossed_x2 = (self.state[:, 0] > (self.bounds[1] - self.size))
        crossed_y1 = (self.state[:, 1] < (self.bounds[2] + self.size))
        crossed_y2 = (self.state[:, 1] > (self.bounds[3] - self.size))
        
        self.state[crossed_x1, 0] = self.bounds[0] + self.size
        self.state[crossed_x2, 0] = self.bounds[1] - self.size

        self.state[crossed_y1, 1] = self.bounds[2] + self.size
        self.state[crossed_y2, 1] = self.bounds[3] - self.size

        self.state[crossed_x1 | crossed_x2, 2] *= -1
        self.state[crossed_y1 | crossed_y2, 3] *= -1

        #chase:
        dist = np.sqrt(((self.state[0,0]-self.state[1,0])**2)+((self.state[0,1]-self.state[1,1])**2))
        self.state[1, :2] += dt*self.state[1, 2:]
        self.state[1, 2] += dt* (self.state[0,0]-self.state[1,0])/dist
        self.state[1,3] += dt* (self.state[0,1]-self.state[1,1])/dist

        #run:
        rand_angle = random.uniform(-1.57, 1.57)
        self.state[0, :2] += dt*self.state[0, 2:]
        self.state[0, 2] += dt*(-1)*((np.cos(rand_angle)*(self.state[0,0] - self.state[1,0]))+(np.sin(rand_angle)*(self.state[0,1]-self.state[1,1])))/dist
        self.state[0, 3] += dt*(-1)*((np.sin(rand_angle)*(self.state[0,1] - self.state[1,1]))+(np.cos(rand_angle)*(self.state[0,0]-self.state[1,0])))/dist
        print(rand_angle)




        

        
    
        
    
# set up initial state

init_state =  ([[0,1,1,0],
                [-3,-3,0,0]])


box = ParticleBox(init_state, size = 0.05)
dt = 1/60.0

global listx,listy
listx=list()
listy=list()

#set up animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect = 'equal', autoscale_on = False,
                     xlim=(-6, 6), ylim=(-6, 6))

fig2 = plt.figure()

fig2.subplots_adjust(left=0, right=1, bottom=0, top=1)

ax2 = fig2.add_subplot(111, aspect = 'equal', autoscale_on = False,
                     xlim=(0, 10), ylim=(-6, 6))


particles, = ax.plot([],[],'bo')
rect = plt.Rectangle(box.bounds[::2],
                     box.bounds[1] - box.bounds[0],
                     box.bounds[3] - box.bounds[2],
                     ec = 'none', lw=2, fc = 'none')
ax.add_patch(rect)
def init():
    global box, rect
    particles.set_data([],[])
    rect.set_edgecolor('none')
    return particles, rect


def animate(i):
    global box, rect, dt, ax, fig
    box.chase(dt)
    listx.append(np.ndarray.tolist(box.state[:, 0])[0])
    listy.append(np.ndarray.tolist(box.state[:, 1])[0])
    ms = int(fig.dpi * 2 * box.size * fig.get_figwidth()/np.diff(ax.get_xbound())[0])
    rect.set_edgecolor('k')
    particles.set_data(box.state[:, 0], box.state[:, 1])    
    particles.set_markersize(ms)
    return particles, rect

"""
x = np.asarray(listx[:1000], dtype = float)
y = np.asarray(listy[:1000], dtype = float)
"""
x = np.asarray([1,2])
y = np.asarray([2,3])
ax2.plot(x,y)
ani = animation.FuncAnimation(fig, animate, frames = 600,
                              interval = 10, blit=True, init_func = init) 

plt.show()     

