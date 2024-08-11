# Ocean Simulation

## Overview

![image](https://github.com/user-attachments/assets/9eaf08d4-44bd-4ace-a177-08328e020193)
Watch the demo video: (https://youtu.be/mi3qrnqoHJg)

Using the water shader, the ocean simulation is played in real time, featured with caustics effect. By clicking on the ocean, ripples can be observe. We also implemented the rotation of the sun and moon, which affects the environment lighting.

## Caustics
﻿﻿﻿![image](https://github.com/user-attachments/assets/ea2bf50c-e38c-4bf6-9d4f-1bc24af8d587)
Caustics is radiosity energy scattered by reflection & refraction.
When an area of light hits the surface, it is scattered or concentrated, and its area projected to the bottom changes. 
We already have the ray tracer, so we can just use it to simulate triangle unit area of light, defined by three points, and thus three rays shot from light source in parallel, and they changes their directions through the water surface, and finally hits on the ground.

![image](https://github.com/user-attachments/assets/cd01d65f-4647-4b1f-b01e-a38bacb175c9)
Parameters:
Area coefficient - definition;
Project length - contrast.




## Workload Distribution

##### Noya: 

- Gerstner wave function

- Circle wave function

- Sun & moon rotation

- Perlin noise for sand

- Correct shadow under transparent surface




##### Tsingtao: 

- Read Noya’s code

- Caustic Effect

- Click interaction

- Circle wave function

- Ripple spread & fade out

- Failed attempt on ray tracing a ton triangles in each shader




