/***************
 * Make Planets*
 * matt h      *
 * ************/

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "maths_funcs.h"
#include "gl_utils.h"
#define _USE_MATH_DEFINES
#define GL_LOG_FILE "gl.log"
#define M_PI 3.14159265359

//window size
int g_gl_width = 1600;
int g_gl_height = 900;
GLFWwindow* g_window = NULL;

struct solarObject {
    GLfloat* points;
    GLint* faces;
    mat4 model_mat;
    //GLuint vao;
} sun, merc, ven, ear, moon, mars, jup, sat, ura, nep ;

int countLabel(std::string mName, char const *label){
    int numLab = 0;
    FILE *stream;
    stream = fopen(mName.c_str(), "r");
    char buf[128];
    while(fscanf(stream, "%s", buf) != EOF){
        if(strcmp(buf,label) == 0)
            numLab++;
    }

    printf("Model has %d %s\n", numLab, label);
    fclose(stream);
    return numLab;
}

void loadVerts(std::string mName, GLfloat verts[]){
    printf("Loading vertices\n");
    int numvert=0;

    FILE *stream;
    stream = fopen(mName.c_str(), "r");

    char buf[128];
    float maxX = -100000000.0f;
    float maxY = -100000000.0f;
    float maxZ = -100000000.0f;
    float minX = 100000000.0f;
    float minY = 100000000.0f;
    float minZ = 100000000.0f;
    float a, b, c;
    char const *label = "v";
    while(fscanf(stream, "%s", buf) != EOF){
        if(strcmp(buf,label) == 0){
            fscanf(stream, "%f %f %f\n", &a, &b, &c);
            //printf("a,b,c: %f, %f, %f", a, b, c);
            
            if(a > maxX){
                maxX = a;
            }
            if(a < minX){
                minX = a;
            }
            if(b > maxY){
                maxY = b;
            }
            if(b < minY){
                minY = b;
            }
            if(c > maxZ){
                maxZ = c;
            }
            if(c < minZ){
                minZ = c;
            }
            
            verts[3*numvert + 0] = 1.0*a;
            verts[3*numvert + 1] = 1.0*b;
            verts[3*numvert + 2] = 1.0*c;
            numvert++;
        }
    }

    float scaleX = maxX-minX;
    float scaleY = maxY-minY;
    float scaleZ = maxZ-minZ;
    float transX = 0.5f*(maxX+minX);
    float transY = 0.5f*(maxY+minY);
    float transZ = 0.5f*(maxZ+minZ);
    printf("Max's: %f, %f, %f\n", maxX, maxY, maxZ);
    printf("Min's: %f, %f, %f\n", minX, minY, minZ);
    printf("Scales: %f, %f, %f\n", scaleX, scaleY, scaleZ);
    printf("Trans: %f, %f, %f\n", transX, transY, transZ);

    for(int i =0;i<numvert;i++){
        verts[3*i+0] = (verts[3*i+0] - transX)/scaleX;
        verts[3*i+1] = (verts[3*i+1] - transY)/scaleY;
        verts[3*i+2] = (verts[3*i+2] - transZ)/scaleZ;
    }
   
    fclose(stream);
    printf("Done loading vertices\n");
}

void loadFaces(std::string mName, GLint faces[]){
    printf("Loading faces...\n");
    int numF = 0;

    FILE *stream;
    stream = fopen(mName.c_str(), "r");

    char buf[128];
    int a,b,c;
    char const *label = "f";
    while(fscanf(stream, "%s", buf) != EOF){
        if(strcmp(buf, label) == 0){
            fscanf(stream, "%d %d %d", &a, &b, &c);
            faces[3*numF+0] = a-1;
            faces[3*numF+1] = b-1;
            faces[3*numF+2] = c-1;
            numF++;
        }
    }
    fclose(stream);
    printf("Done loading faces\n");
}

void computeFaceNormals(GLfloat faceNormals[], GLfloat verts[], GLint faces[], int numFaces){
	for (int i = 0; i < numFaces; i++){
		int idx1 = faces[i*3 + 0];
		int idx2 = faces[i*3 + 1];
		int idx3 = faces[i*3 + 2];
		float ux = verts[idx1*3 + 0] - verts[idx2*3 + 0];
		float uy = verts[idx1*3 + 1] - verts[idx2*3 + 1];
		float uz = verts[idx1*3 + 2] - verts[idx2*3 + 2];
		float vx = verts[idx1*3 + 0] - verts[idx3*3 + 0];
		float vy = verts[idx1*3 + 1] - verts[idx3*3 + 1];
		float vz = verts[idx1*3 + 2] - verts[idx3*3 + 2];
		float nx = uy*vz - uz*vy;
		float ny = uz*vx - ux*vz;
		float nz = ux*vy - uy*vx;
		float mag = sqrt(nx*nx + ny*ny + nz*nz);
		
	//	cout << "avg norm" << nx << ", " << ny << ", " << nz << endl;
		//cout << "mag: " << mag << endl;
		faceNormals[3*i + 0] = nx/mag;
		faceNormals[3*i + 1] = ny/mag;
		faceNormals[3*i + 2] = nz/mag;
	}
	
}

void computeVertNormals(GLfloat normals[], GLfloat verts[], int numVerts, GLint faces[], int numFaces, GLfloat faceNormals[]){
	for (int i = 0; i < numVerts; i++){
		float avgX = 0.0;
		float avgY = 0.0;
		float avgZ = 0.0;
		int numF_vert = 0;
		
		//find all the faces that contain this vertex
		for (int j = 0; j < numFaces; j++){
			int found = 0;
			for (int k = 0; k < 3; k++)
				if (faces[j*3 + k] == i)
					found = 1;
			if (found){
				avgX += faceNormals[j*3 + 0];
				avgY += faceNormals[j*3 + 1];
				avgZ += faceNormals[j*3 + 2];				
				numF_vert++;				
			}			
		}
		
		avgX /= numF_vert;
		avgY /= numF_vert;
		avgZ /= numF_vert;
		//cout << "avg norm" << avgX << ", " << avgY << ", " << avgZ << endl;
		normals[i*3 + 0] = avgX;
		normals[i*3 + 1] = avgY;
		normals[i*3 + 2] = avgZ;
	}
}

int main(){
    
    restart_gl_log();
    start_gl();

    glewExperimental = GL_TRUE;
    glewInit();

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    /** function calls for creating the planet goes here **/
    std::string modelName = 
		"C:/Users/Gabe/Documents/Visual Studio 2015/Projects/Gallagher_370_Final/Gallagher_370_Final/MattProject/MattProject/sphere.obj";
    
    int v = countLabel(modelName, "v");
    GLfloat* verts = new GLfloat[3*v];
    loadVerts(modelName, verts);
	
    int numFaces = countLabel(modelName, "f");
    GLint* faces = new GLint[3*numFaces];
    loadFaces(modelName, faces);
	
    GLfloat* faceNormals = new GLfloat[3*numFaces];
    computeFaceNormals(faceNormals, verts, faces, numFaces);
	
    GLfloat* vertNormals = new GLfloat[3*v];
    computeVertNormals(vertNormals, verts, v, faces, numFaces, faceNormals);

    GLfloat* points = new GLfloat[9*numFaces];
    GLfloat* normals = new GLfloat[9*numFaces];
	for (int i = 0; i < numFaces; i++){
        int idx1 = faces[3*i + 0];
        int idx2 = faces[3*i + 1];
        int idx3 = faces[3*i + 2];
        points[i*9 + 0] = verts[3*idx1+0];
        points[i*9 + 1] = verts[3*idx1+1];
        points[i*9 + 2] = verts[3*idx1+2];
        points[i*9 + 3] = verts[3*idx2+0];
        points[i*9 + 4] = verts[3*idx2+1];
        points[i*9 + 5] = verts[3*idx2+2];
        points[i*9 + 6] = verts[3*idx3+0];
        points[i*9 + 7] = verts[3*idx3+1];
        points[i*9 + 8] = verts[3*idx3+2];
        normals[i*9 + 0] = vertNormals[3*idx1+0];
        normals[i*9 + 1] = vertNormals[3*idx1+1];
        normals[i*9 + 2] = vertNormals[3*idx1+2];
        normals[i*9 + 3] = vertNormals[3*idx2+0];
        normals[i*9 + 4] = vertNormals[3*idx2+1];
        normals[i*9 + 5] = vertNormals[3*idx2+2];
        normals[i*9 + 6] = vertNormals[3*idx3+0];
        normals[i*9 + 7] = vertNormals[3*idx3+1];
        normals[i*9 + 8] = vertNormals[3*idx3+2];
        }
	int numPoints = 3*numFaces; 

    sun.points = points;
    sun.faces = faces;
    merc.points = points;
    merc.faces = faces;
    
    GLuint points_vbo;
    glGenBuffers(1, &points_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * numPoints * sizeof(GLfloat), points, GL_STATIC_DRAW);

    GLuint normals_vbo;
    glGenBuffers(1, &normals_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
    glBufferData(GL_ARRAY_BUFFER, 9 * numFaces * sizeof(GLfloat), normals, GL_STATIC_DRAW);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    //sun.vao = vao;
    //merc.vao = vao;

    GLuint shader = create_programme_from_files(
		"C:/Users/Gabe/Documents/Visual Studio 2015/Projects/Gallagher_370_Final/Gallagher_370_Final/MattProject/MattProject/vs.glsl",
		"C:/Users/Gabe/Documents/Visual Studio 2015/Projects/Gallagher_370_Final/Gallagher_370_Final/MattProject/MattProject/fs.glsl"
	);

    GLuint mercShader = create_programme_from_files(
		"C:/Users/Gabe/Documents/Visual Studio 2015/Projects/Gallagher_370_Final/Gallagher_370_Final/MattProject/MattProject/merc_vs.glsl",
		"C:/Users/Gabe/Documents/Visual Studio 2015/Projects/Gallagher_370_Final/Gallagher_370_Final/MattProject/MattProject/merc_fs.glsl"
	);

    //stolen from anton
    
    // input variables
    float near = 0.1f; //clipping plane
    float far = 200.0f; //clipping plane
    float fov = 67.0f * ONE_DEG_IN_RAD; // 67 degrees to radians
    float aspect = (float)g_gl_width / (float)g_gl_height;

    // matrix components
    float range = tan(fov*0.5f) * near;
    float Sx = (2.0f * near) / (range * aspect + range * aspect);
    float Sy = near / range;
    float Sz = -(far + near) / (far - near);
    float Pz = -(2.0f * far * near) / (far - near);

    GLfloat proj_mat[] = {
        Sx, 0.0f, 0.0f, 0.0f,
        0.0f, Sy, 0.0f, 0.0f,
        0.0f, 0.0f, Sz, -1.0f,
        0.0f, 0.0f, Pz, 0.0f
    };

    // create view matrix
    float cam_pos[] = {0.0f, 0.0f, 2.0f};
    float cam_yaw = 0.0f;
    mat4 T = translate(identity_mat4(), vec3 (-cam_pos[0], -cam_pos[1], -cam_pos[2]));
    mat4 R = rotate_y_deg(identity_mat4(), -cam_yaw);
    mat4 view_mat = R*T;

    //mat4 model_mat = identity_mat4();

    glUseProgram(shader);
    int view_mat_location = glGetUniformLocation(shader, "view_mat");
    glUniformMatrix4fv(view_mat_location, 1, GL_FALSE, view_mat.m);
    int proj_mat_location = glGetUniformLocation(shader, "projection_mat");
    glUniformMatrix4fv(proj_mat_location, 1, GL_FALSE, proj_mat);
    //int model_mat_location = glGetUniformLocation(shader, "model_mat");
    //glUniformMatrix4fv(model_mat_location, 1, GL_FALSE, model_mat.m);
    int sun_mat_location = glGetUniformLocation(shader, "sun_model_mat");
    glUniformMatrix4fv(sun_mat_location, 1, GL_FALSE, sun.model_mat.m);

    glUseProgram(mercShader);
    int merc_view_mat_location = glGetUniformLocation(mercShader, "view_mat");
    glUniformMatrix4fv(merc_view_mat_location, 1, GL_FALSE, view_mat.m);
    int merc_proj_mat_location = glGetUniformLocation(mercShader, "projection_mat");
    glUniformMatrix4fv(merc_proj_mat_location, 1, GL_FALSE, proj_mat);
    int merc_mat_location = glGetUniformLocation(mercShader, "model_mat");
    glUniformMatrix4fv(merc_mat_location, 1, GL_FALSE, merc.model_mat.m);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);


    while(!glfwWindowShouldClose(g_window)){

        _update_fps_counter (g_window);
        double current_seconds = glfwGetTime();

        sun.model_mat = rotate_y_deg(identity_mat4(), current_seconds*50.0f);
        merc.model_mat = rotate_y_deg(identity_mat4(), current_seconds*15.0f) * 
            translate(identity_mat4(), vec3(-1.5, 0.0, 0.0)) *
            scale(identity_mat4(), vec3(0.2, 0.2, 0.2));

        //and then the rest of the solar system objects
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, g_gl_width, g_gl_height);
        glUseProgram(shader);
        //glUniformMatrix4fv(model_mat_location, 1, GL_FALSE, model_mat.m);
        glUniformMatrix4fv(sun_mat_location, 1, GL_FALSE, sun.model_mat.m);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, numPoints);
        glUseProgram(mercShader);
        glUniformMatrix4fv(merc_mat_location, 1, GL_FALSE, merc.model_mat.m);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, numPoints);
        glfwPollEvents();
        if(GLFW_PRESS == glfwGetKey(g_window, GLFW_KEY_ESCAPE)){
            glfwSetWindowShouldClose(g_window, 1);
        }
        glfwSwapBuffers(g_window);
    }
    glfwTerminate();
    return 0;
}