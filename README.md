# GPU-Raytracer
This is a GPU based Raytracer using C++/GLSL. ALL object defining and ray tracing take place in shaders. To render basic objects, use color shader. To render ocean simulation, use water shader. No additional libraries needed.

## Raytracer

Ray-tracing is done in the vertex and fragment shaders: GI.vs and GI.fs. The majority of the tracing is done in the fragment shader and the vertex shader is only passing the value of screen coordinates to fragment shader. 

![Final](Image/Raytracing.bmp)

## Tone Reproduction

Tone Reproduction is done in the main.cpp in the src folder in function tone_reproduction() and tone_reproduction_Advanced(). Uncomment line 376 and 377 to enable both tone reproduction functions. Make sure to close the picture window before closing the console because tone reproduction only runs after picture window is closed.

## Ocean Simulation
