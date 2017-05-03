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
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "obj_model.h"
#include <vector>

using namespace std;

//window size
int g_gl_width = 1600;
int g_gl_height = 900;
GLFWwindow* g_window = NULL;

bool load_texture(const char *file_name, GLuint *tex ) {
	int x, y, n;
	int force_channels = 4;
	unsigned char *image_data = stbi_load( file_name, &x, &y, &n, force_channels );
	if ( !image_data ) {
		fprintf( stderr, "ERROR: could not load %s\n", file_name );
		return false;
	}
	// NPOT check
	if ( ( x & ( x - 1 ) ) != 0 || ( y & ( y - 1 ) ) != 0 ) {
		fprintf( stderr, "WARNING: texture %s is not power-of-2 dimensions\n",
						 file_name );
	}

	int width_in_bytes = x * 4;
	unsigned char *top = NULL;
	unsigned char *bottom = NULL;
	unsigned char temp = 0;
	int half_height = y / 2;

	for ( int row = 0; row < half_height; row++ ) {
		top = image_data + row * width_in_bytes;
		bottom = image_data + ( y - row - 1 ) * width_in_bytes;
		for ( int col = 0; col < width_in_bytes; col++ ) {
			temp = *top;
			*top = *bottom;
			*bottom = temp;
			top++;
			bottom++;
		}
	}

	glGenTextures(1, tex);
	glActiveTexture(0);
	glBindTexture(GL_TEXTURE_2D, *tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, x, y, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);

	glGenerateMipmap( GL_TEXTURE_2D );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
	GLfloat max_aniso = 0.0f;
	glGetFloatv( GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &max_aniso );
	// set the maximum!  glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, max_aniso ); 
	return true;
}

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

void load_texture_coords(string fileName, GLfloat* texs) {
  printf("Loading texture coordinates\n");
  int numCoords = 0;

  FILE* stream;
  stream = fopen(fileName.c_str(), "r");

  char buf[128];
  char const* label = "vt";
  float a, b;
  float maxX = -100000000.0f;
  float maxY = -100000000.0f;
  float minX = 100000000.0f;
  float minY = 100000000.0f;

  while (fscanf(stream, "%s", buf) != EOF) {
    if (strcmp(buf, label) == 0) {
      fscanf(stream, "%f %f\n", &a, &b);
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
            
      texs[2 * numCoords] = a * 1.0;
      texs[2 * numCoords + 1] = b * 1.0;
      numCoords++;
    }
  }

  fclose(stream);
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
            fscanf(stream, "%d%*s %d%*s %d%*s", &a, &b, &c);
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
  cout << "Computing vert norms" << endl;
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
    std::string modelName = "tex_sphere.obj";
    
    int v = countLabel(modelName, "v");
    GLfloat* verts = new GLfloat[3 * v];
    loadVerts(modelName, verts);
	
    int numFaces = countLabel(modelName, "f");
    GLint* faces = new GLint[3 * numFaces];
    loadFaces(modelName, faces);
	
    GLfloat* faceNormals = new GLfloat[3 * numFaces];
    computeFaceNormals(faceNormals, verts, faces, numFaces);
	
    GLfloat* vertNormals = new GLfloat[3 * v];
    computeVertNormals(vertNormals, verts, v, faces, numFaces, faceNormals);

    GLfloat* points = new GLfloat[9 * numFaces];
    GLfloat* normals = new GLfloat[9 * numFaces];

	  for (int i = 0; i < numFaces; i++) {
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

    // Create the pointers for the celestial objects
    ObjModel* sun_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* merc_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* venus_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* earth_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* mars_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* jupiter_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* saturn_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* uranus_obj = new ObjModel(points, faces, NULL, numPoints, numFaces, 0);
    ObjModel* neptune_obj =  new ObjModel(points, faces, NULL, numPoints, numFaces, 0);

    // Create a vector and add the celestial objects
    vector<ObjModel*> celestials = {};
    celestials.push_back(sun_obj);
    celestials.push_back(merc_obj);
    celestials.push_back(venus_obj);
    celestials.push_back(earth_obj);
    celestials.push_back(mars_obj);
    celestials.push_back(jupiter_obj);
    celestials.push_back(saturn_obj);
    celestials.push_back(uranus_obj);
    celestials.push_back(neptune_obj);

    int coord_count = countLabel(modelName, "vt");
    GLfloat* tex_coords = new GLfloat[coord_count * 2];
    load_texture_coords(modelName, tex_coords);

    GLuint points_vbo;
    glGenBuffers(1, &points_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * numPoints * sizeof(GLfloat), points, GL_STATIC_DRAW);

    GLuint normals_vbo;
    glGenBuffers(1, &normals_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * numPoints * sizeof(GLfloat), normals, GL_STATIC_DRAW);

    GLuint tex_vbo;
    glGenBuffers(1, &tex_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, tex_vbo);
    glBufferData(GL_ARRAY_BUFFER, 2 * numPoints * sizeof(GLfloat), tex_coords, GL_STATIC_DRAW);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    glBindBuffer(GL_ARRAY_BUFFER, tex_vbo);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    sun_obj->set_shader(create_programme_from_files("vs.glsl", "fs.glsl"));
    merc_obj->set_shader(create_programme_from_files("merc_vs.glsl", "merc_fs.glsl"));
    venus_obj->set_shader(create_programme_from_files("merc_vs.glsl", "venus_fs.glsl"));
    earth_obj->set_shader(create_programme_from_files("merc_vs.glsl", "earth_fs.glsl"));
    mars_obj->set_shader(create_programme_from_files("merc_vs.glsl", "mars_fs.glsl"));
    jupiter_obj->set_shader(create_programme_from_files("merc_vs.glsl", "jupiter_fs.glsl"));
    saturn_obj->set_shader(create_programme_from_files("merc_vs.glsl", "saturn_fs.glsl"));
    uranus_obj->set_shader(create_programme_from_files("merc_vs.glsl", "uranus_fs.glsl"));
    neptune_obj->set_shader(create_programme_from_files("merc_vs.glsl", "neptune_fs.glsl"));

    //stolen from anton
    
    // input variables
    float near = 0.1f; //clipping plane
    float far = 500.0f; //clipping plane
    float fov = 80.0f * ONE_DEG_IN_RAD; // 67 degrees to radians
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
    float cam_pos[] = {0.0f, 2.0f, 4.0f};
    float cam_yaw = -35.0f;

    mat4 T = translate(identity_mat4(), vec3 (-cam_pos[0], -cam_pos[1], -cam_pos[2]));
    mat4 R = rotate_x_deg(identity_mat4(), -cam_yaw);
    mat4 view_mat = R*T;

    for (ObjModel* model : celestials) {
      glUseProgram(model->get_shader());

      int view_mat_location = glGetUniformLocation(model->get_shader(), "view_mat");
      glUniformMatrix4fv(view_mat_location, 1, GL_FALSE, view_mat.m);

      int proj_mat_location = glGetUniformLocation(model->get_shader(), "projection_mat");
      glUniformMatrix4fv(proj_mat_location, 1, GL_FALSE, proj_mat);

      int mat_location = glGetUniformLocation(model->get_shader(), "model_mat");
      glUniformMatrix4fv(mat_location, 1, GL_FALSE, model->get_model_mat().m);

      model->set_mat_location(mat_location);
    }

    GLuint tex;
    (load_texture( "mercurymap.jpg", &tex));

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);

    double earth_size = 0.075;
    double earth_distance = -3.5;
    float earth_speed = 30.0f;

    // Celestials configuration
    sun_obj->set_rot_speed(0.0, 1, 0.0);
    sun_obj->set_translation(0.0, 0.0, 0.0);

    merc_obj->set_scale(earth_size * 0.383, earth_size * 0.383, earth_size * 0.383);
    merc_obj->set_rot_speed(0.0, 1.59 * earth_speed, 0.0);
    merc_obj->set_translation(earth_distance * 0.387, 0.0, 0.0);

    venus_obj->set_scale(earth_size * 0.949, earth_size * 0.949, earth_size * 0.949);
    venus_obj->set_rot_speed(0.0, 1.18 * earth_speed, 0.0);
    venus_obj->set_translation(earth_distance * 0.723, 0.0, 0.0);

    earth_obj->set_scale(earth_size, earth_size, earth_size);
    earth_obj->set_rot_speed(0.0, earth_speed, 0.0);
    earth_obj->set_translation(earth_distance, 0.0, 0.0);

    mars_obj->set_scale(earth_size * 0.532, earth_size * 0.532, earth_size * 0.532);
    mars_obj->set_rot_speed(0.0, 0.808 * earth_speed, 0.0);
    mars_obj->set_translation(earth_distance * 1.52, 0.0, 0.0);

    jupiter_obj->set_scale(earth_size * 11.21, earth_size * 11.21, earth_size * 11.21);
    jupiter_obj->set_rot_speed(0.0, 0.439 * earth_speed, 0.0);
    jupiter_obj->set_translation(earth_distance * 5.2, 0.0, 0.0);

    saturn_obj->set_scale(earth_size * 9.45, earth_size * 9.45, earth_size * 9.45);
    saturn_obj->set_rot_speed(0.0, 0.325 * earth_speed, 0.0);
    saturn_obj->set_translation(earth_distance * 9.58, 0.0, 0.0);
    
    uranus_obj->set_scale(earth_size * 4.01, earth_size * 4.01, earth_size * 4.01);
    uranus_obj->set_rot_speed(0.0, 0.228 * earth_speed, 0.0);
    uranus_obj->set_translation(earth_distance * 19.20, 0.0, 0.0);

    neptune_obj->set_scale(earth_size * 3.88, earth_size * 3.88, earth_size * 3.88);
    neptune_obj->set_rot_speed(0.0, 0.182 * earth_speed, 0.0);
    neptune_obj->set_translation(earth_distance * 30.05, 0.0, 0.0);

    while(!glfwWindowShouldClose(g_window)){

        _update_fps_counter (g_window);
        double current_seconds = glfwGetTime();

        // Iterate through and update all the celestial bodies.
        for (ObjModel* model : celestials) {
          model->update(current_seconds);
        }
       
        //and then the rest of the solar system objects
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, g_gl_width, g_gl_height);

        for (ObjModel* model : celestials) {
          // Swap to the model's shader
          glUseProgram(model->get_shader());
          glUniformMatrix4fv(model->get_mat_location(), 1, GL_FALSE, model->get_model_mat().m);
          glBindVertexArray(vao);
          glDrawArrays(GL_TRIANGLES, 0, numPoints);
        }

        glfwPollEvents();

        if(GLFW_PRESS == glfwGetKey(g_window, GLFW_KEY_ESCAPE)){
            glfwSetWindowShouldClose(g_window, 1);
        }

        glfwSwapBuffers(g_window);
    }
    glfwTerminate();
    return 0;
}

