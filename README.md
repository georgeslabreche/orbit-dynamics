R0008R - Introduction to Space Mechanics and Electronics
===- Teacher: Leonard Felicetti.
- Assignment 2: Orbit Dynamics.


Task 1
------
After having written a matlab code that performs the integration of the equation of motion of the orbit dynamics of a spacecraft around Earth, analyze and illustrate the orbital motion obtained with the following initial orbital parameters:

| Name | LEO  |  GTO  |   GEO   | Molniya | Tundra  |  MEO  | 
|:----:|-----:|------:|--------:|--------:|--------:|------:|
|  𝑎   | 7978 | 24400 |Calculate|  26555  |  42164  | 26560 |
|  𝑒   | 0.05 | 0.73  |   0.2   |   0.72  |  0.075  | 0.01  |
|  𝑖   | 47.2 |   7   |   10    |   63.4  |    43   |  53.7 |
|  Ω   | 30  |   10  |   30    |   180   |    195  |  145.5 |
|  𝜔   |  45 |   5   |   325   |   270   |  270    |   52.4 |
|  𝑀  |70   |  250  |   300   |   90    |   305   |  -37.3 |

Set the simulation in such way that at time 𝑡<sub>0</sub>=0 the position of the satellite in the orbit is given by the mean anomaly 𝑀 in the last row of the table. Calculate the initial position and velocity of the satellite and integrate its motion for 10 orbital periods. Plot the trajectory of the spacecraft around the Earth and the behavior of the orbital parameters during the simulation time. Provide comments on all the figures and give a physical interpretation even in relation of the kind of orbit that you are analyzing.

Task 2
------
Plot the satellite ground tracks of the satellites analyzed in the previous task. Assume that the Greenwich meridian (longitude 𝜆=0 𝑑𝑒𝑔) is aligned with the first axis (𝛾 point of Aries) of the ECI reference frame at 𝑡0=0. Describe and give a physical interpretation of the obtained plots even in relation on the kind of orbit you are analyzing.

Task 3
------
Four ground stations are located in different positions around the Earth, as shown in the following table:

|Ground Station|   Kiruna  |  Malindi  |Cape Canaveral| Tanegashima| 
|:------------:|----------:|----------:|-------------:|-----------:|
|  Latitude    |67°53'22"𝑁 |2°56'18"𝑁 |28°31'26"𝑁    |30°23'60"𝑁 |
|  Longitude   |20°13'06"𝐸 |40°12'45"𝐸 |80°39'3"𝑊    |130°58'7"𝐸  |

For each of the satellites analyzed on the previous cases, check the visibility condition and plot the elevation and azimuth angles for which the satellite can be seen from these stations during the simulation period. Describe and give an interpretation of the plots you obtain.

Task 4
------
A satellite is initially orbiting in equatorial LEO orbit at 300𝑘𝑚 altitude, 𝑒=0.001, 𝑖=0.01 𝑑𝑒𝑔, Ω=10 𝑑𝑒𝑔, 𝜔=270 𝑑𝑒𝑔. After 2 orbits, the satellites performs an orbital transfer from initial orbit (LEO) to a geostationary orbit (GEO) and then it stays in GEO orbit for other 2 orbital periods.

**a)** After having calculated the theoretical values of Δ𝑉 to be provided at perigee and apogee of the transfer orbit, please simulate the scenario in order to obtain the plot of the orbital parameters

**b)** By assuming that the thrust is provided impulsively by a bipropellant thruster with a specific impulse 𝐼𝑠𝑝=300𝑠 and the initial mass of the spacecraft is 𝑚<sub>0</sub>=2500𝑘𝑔, calculate the propellant mass necessary to complete the transfer.**c)** If the initial orbit has an inclination 𝑖=10 𝑑𝑒𝑔 (the other initial parameters do not change), the orbit transfer has to include an orbital plane change maneuver. Select the point of the transfer orbit that is preferable for performing such orbit plane change maneuver, saving propellant mass. Explain the reason of your choice by analyzing the problem quantitatively through calculations and simulations. Simulate the new scenario in order to obtain the plots of the orbital parameters and visualization the three obtained orbits. Comment and interpret the obtained results.