#version 330 core
layout (location = 0) in vec3 aPos;

layout (location = 1) in vec3 aColor;

out vec3 ourColor;
out vec2 screenCoord;
//out float time_t;
//uniform float time_t;

void main()
{
    gl_Position = vec4(aPos, 1.0f); 
    //gl_Position = vec4(aPos, 1.0f);
    ourColor = aColor;
    screenCoord = (vec2(aPos.x, aPos.y) + 1.0) / 2.0;
}