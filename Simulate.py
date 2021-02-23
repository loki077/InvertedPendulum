import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

#Define constants that go into simulation
massCart = 5.0; #[kg]
massPendulum = 1.0; #[kg]
Length = 0.5; # [m]
g = 9.82; #[m/s^2]
Force = 0; #[N]
Torque = 0;
Bcart = 2; #[N / m/s]
BPendulum = 0.1; #[Nm / rad/s]

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


	Theta2Dot = ( Torque_eq/(massPendulum*Length)*(massPendulum+massCart) - Force_eq*math.cos(InitialValues[0]) + g*math.sin(InitialValues[0])*(massCart+massPendulum) \
		-massPendulum*Length*math.sin(InitialValues[0])*math.cos(InitialValues[0])*InitialValues[1]**2)/( Length*(massCart+massPendulum*math.sin(InitialValues[0])) );
	x2Dot = ( -massPendulum*math.sin(InitialValues[0])*(g*math.cos(InitialValues[0]) -Length*InitialValues[1]**2) -  \
		Torque_eq/Length*math.cos(InitialValues[0])+Force_eq  ) / ( massCart+massPendulum*math.sin(InitialValues[0]) );

	#Theta2Dot = ( (massCart+massPendulum) * g * math.sin(InitialValues[0]) - massPendulum*Length*math.sin(InitialValues[0]) * math.cos(InitialValues[0]) * InitialValues[1]**2 \
	#	-(Force-Bcart*InitialValues[3])*math.cos(InitialValues[0]) ) / (Length * (massCart+massPendulum*(math.sin(InitialValues[0]))**2));
	#x2Dot = ((Force-Bcart*InitialValues[3]) + massPendulum*math.sin(InitialValues[0]) * (Length*InitialValues[1]-g * math.cos(InitialValues[0])))/(massCart+massPendulum*(math.sin(InitialValues[0]))**2)


	return Theta1Dot, Theta2Dot, x1Dot, x2Dot

SimValues = integrate.odeint(PendulumOnCartSim, InitialValues, Time, args=Args)

#plot figures
plt.plot(Time,SimValues[:,3])
plt.show()