#ifndef FINAL_OBJMODEL_H
#define FINAL_OBJMODEL_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>
#include <vector>

#include "maths_funcs.h"

/**
 * Helper class that represents an .obj model
 */
class ObjModel {
  private:
    // Global variables
    vec3                 rot_speed;
    vec3                 rot;
    vec3                 obj_scale;
    vec3                 translation;
    std::vector<GLfloat> points;
    std::vector<GLint>   faces;
    std::vector<GLfloat> texture_coords;
    mat4                 model_mat;
    GLuint               shader;
    int                  mat_location;

    template<typename T> std::vector<T> convert_to_vector(T*, int);
    template<typename T> T* convert_to_array(std::vector<T>);
  public:
    // Default constructor
    ObjModel(GLfloat* points, GLint* faces, GLfloat* tex_coords,
        int point_count, int face_count, int tex_count);

    // Void methods
    void update(double seconds);

    // Getter methods
    GLfloat*   get_points();
    GLint*     get_faces();
    GLfloat*   get_texture_coords();
    mat4       get_model_mat();
    GLuint     get_shader();
    int        get_mat_location();

    // Setter methods
    void set_rot(float r_x, float r_y, float r_z);
    void set_rot_speed(float speed_x, float speed_y, float speed_z);
    void set_scale(double s_x, double s_y, double s_z);

    void set_translation(float x, float y, float z);

    void set_shader(GLuint shader);

    void set_mat_location(int);
};

#endif
