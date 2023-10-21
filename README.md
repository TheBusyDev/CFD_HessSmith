# Hess Smith Panel Method
A simple *CFD* method to compute __ideal flows__ around bodies.\
The core idea is to approximate the boundary of the body as a series of *sinks*, *sources* and *vortices*, computing the volume flow rate or the circulation for each of them. To complete these calculations, the *Kutta condition* is prescribed for the trailing edge.\
Then it's possible to find the __velocity field__ around the body. Here is an interesting example of the flow surrounding a __NACA arifoil__ and a __cylinder__ (obviously, since we're considering an ideal flow, the wake is NOT simulated, so we cannot rely on this model to estimate the behaviour at the rear of the body).\
\
![](/images/airfoil_u.png)
![](/images/airfoil_v.png)
\
\
![](/images/cylinder_u.png)
![](/images/cylinder_v.png)
