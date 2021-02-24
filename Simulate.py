import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

#Define constants that go into simulation
massCart = 5; #[kg]
massPendulum = 1.0; #[kg]
Length = 1; # [m]
g = 9.82; #[m/s^2]
Force = 0; #[N]
Torque = 0;
Bcart = 1; #[N / m/s]
BPendulum = 0.4; #[Nm / rad/s]

Args =  (massCart, Length, g, massPendulum, Bcart, Force, Torque);

Ts = 0.01; #[s]
SimulationTime = 20.0; #[s]

InitialAngle = math.pi/2; #[rad]
InitialAngleRate = 0.0; #[rad/s]
InitialPosition = 0.0; #[m]
InitialVelocity = 0.0; #[m]
InitialValues = (InitialAngle, InitialAngleRate, InitialPosition, InitialVelocity);
#print(InitialValues[3])

#Create Time vector
Time = np.linspace(0, SimulationTime, int(SimulationTime/Ts));
#print(Time)

def PendulumOnCartSim(InitialValues, t, massCart, Length, g, massPendulum, Bcart, Force, Torque):

	Theta1Dot = InitialValues[1];
	x1Dot = InitialValues[3];

	Force_eq = Force-Bcart*InitialValues[3]; #total force acting on system
	Torque_eq = Torque-BPendulum*InitialValues[1]; #total torque acting on system


	Theta2Dot = ( Torque_eq/(massPendulum*Length)*(massPendulum+massCart) - Force_eq*np.cos(InitialValues[0]) + g*np.sin(InitialValues[0])*(massCart+massPendulum) \
		-massPendulum*Length*np.sin(InitialValues[0])*np.cos(InitialValues[0])*InitialValues[1]**2)  /  ( Length*(massCart+massPendulum*np.sin(InitialValues[0])**2) );

	x2Dot = 1/(massCart+massPendulum) * (Force_eq + massPendulum*Length*np.sin(InitialValues[0])*InitialValues[1]**2 - ( Torque_eq/Length*np.cos(InitialValues[0])*(massCart+massPendulum) \
		- Force_eq*massPendulum*np.cos(InitialValues[0])**2 + massPendulum*g*np.sin(InitialValues[0])*np.cos(InitialValues[0])*(massPendulum+massCart) - massPendulum**2*Length*np.sin(InitialValues[0])*np.cos(InitialValues[0])**2*InitialValues[1]**2) \
	/( massCart+massPendulum*np.sin(InitialValues[0])**2 ) )


	return Theta1Dot, Theta2Dot, x1Dot, x2Dot

SimValues = integrate.odeint(PendulumOnCartSim, InitialValues, Time, args=Args)

#Position of cart and pendulum
xPositionCart = SimValues[:,2];
yPositionCart = np.zeros(len(SimValues[:,1]));
#PositionCart = (xPositionCart, yPositionCart);

xPositionPendulum = xPositionCart + Length*np.sin(SimValues[:,0]);
yPositionPendulum = Length*np.cos(SimValues[:,0]);

##plot figures
#plt.plot(xPositionPendulum,yPositionPendulum)
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-3, 3), ylim=(-3, 3))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
	thisx = [xPositionCart[i], xPositionPendulum[i]]
	thisy = [yPositionCart[i], yPositionPendulum[i]]

	line.set_data(thisx, thisy)
	time_text.set_text(time_template % (i*Ts))
	return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(SimValues)),
                              interval=10, blit=True, init_func=init)

plt.show()