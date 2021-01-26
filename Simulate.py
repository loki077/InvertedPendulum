import numpy as np
from scipy import integrate
import matplotlib.pyplot 

#Define constants that go into simulation
massCart = 1.0; #[kg]
massPendulum = 1.0; #[kg]
Length = 1.0; # [m]
g = 9.82; #[m/s^2]
Ts = 0.01; #[s]
SimulationTime = 10.0; #[s]

InitialAngle = 0.0; #[rad]
InitialAngleRate = 0.0; #[rad/s]
InitialPosition = 0.0; #[m]
InitialVelocity = 0.0; #[m]
InitialValues = (InitialAngle, InitialAngleRate, InitialPosition, InitialVelocity);

#Create Time vector
Time = np.linspace(0, SimulationTime, int(SimulationTime/Ts));

def PendulumOnCartSim(InitialValues, massCart, Length, g, Ts, massPendulum, Force):
	Theta1Dot = InitialValues[1];
	x1Dot = InitialValues[3];

	Theta2Dot = ( (massCart+massPendulum) * g * sin(InitialValues[0]) - massPendulum*Length*sin(InitialValues[0]) * cos(InitialValues[0]) * InitialValues[1]**2 \
		-Force*cos(InitialValues[0]) ) / (Length * (massCart+massPendulum*(sin(InitialValues[0]))**2));
	x2Dot = (F + massPendulum*sin(InitialValues[0]) * (Length*InitialValues[1]-g * cos(InitialValues[0])))/(massCart+massPendulum*(sin(InitialValues[0]))**2)


	return Theta1Dot, Theta2Dot, x1Dot, x2Dot

Args =  (InitialValues, Time, massCart, Length, g, Ts, massPendulum, 0);
SimValues = integrate.odeint(PendulumOnCartSim, Time, Args)