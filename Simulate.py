import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import math

#Define constants that go into simulation
massCart = 1.0; #[kg]
massPendulum = 1.0; #[kg]
Length = 1.0; # [m]
g = 9.82; #[m/s^2]
Force = 0;

Args =  (massCart, Length, g, massPendulum, Force);

Ts = 0.01; #[s]
SimulationTime = 10.0; #[s]

InitialAngle = math.pi/2; #[rad]
InitialAngleRate = 0.0; #[rad/s]
InitialPosition = 0.0; #[m]
InitialVelocity = 0.0; #[m]
InitialValues = (InitialAngle, InitialAngleRate, InitialPosition, InitialVelocity);
#print(InitialValues[3])

#Create Time vector
Time = np.linspace(0, SimulationTime, int(SimulationTime/Ts));
#print(Time)

def PendulumOnCartSim(InitialValues, t, massCart, Length, g, massPendulum, Force):

	Theta1Dot = InitialValues[1];
	x1Dot = InitialValues[3];

	Theta2Dot = ( (massCart+massPendulum) * g * math.sin(InitialValues[0]) - massPendulum*Length*math.sin(InitialValues[0]) * math.cos(InitialValues[0]) * InitialValues[1]**2 \
		-Force*math.cos(InitialValues[0]) ) / (Length * (massCart+massPendulum*(math.sin(InitialValues[0]))**2));
	x2Dot = (Force + massPendulum*math.sin(InitialValues[0]) * (Length*InitialValues[1]-g * math.cos(InitialValues[0])))/(massCart+massPendulum*(math.sin(InitialValues[0]))**2)


	return Theta1Dot, Theta2Dot, x1Dot, x2Dot

SimValues = integrate.odeint(PendulumOnCartSim, InitialValues, Time, args=Args)
plt.plot(Time,SimValues)