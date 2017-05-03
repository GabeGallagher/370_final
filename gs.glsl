//GEOMETRY SHADER

#version 450

//A Geometry Shader (GS) is a Shader program written in GLSL that governs the processing 
//of Primitives. Geometry shaders reside between the Vertex Shaders (or the optional 
//Tessellation stage) and the fixed-function Vertex Post-Processing stage.

//The other feature was GS instancing, which allows multiple invocations to operate over 
//the same input primitive. This makes layered rendering easier to implement and possibly 
//faster performing, as each layer's primitive(s) can be computed by a separate GS instance.

//Note: While geometry shaders have had previous extensions like GL_EXT_geometry_shader4 
//and GL_ARB_geometry_shader4, these extensions expose the API and GLSL functionality in 
//very different ways from the core feature. This page describes only the core feature.

layout (triangles) in;
layout (triangle_strip, max_vertices=18) out;

uniform mat4 shadowMatrices[6];

out vec4 FragPos; // FragPos from GS (output per emitvertex)

void main()
{
    for(int face = 0; face < 6; ++face)
    {
        gl_Layer = face; // built-in variable that specifies to which face we render.
        for(int i = 0; i < 3; ++i) // for each triangle's vertices
        {
            FragPos = gl_in[i].gl_Position;
            gl_Position = shadowMatrices[face] * FragPos;
            EmitVertex();
        }    
        EndPrimitive();
    }
} 