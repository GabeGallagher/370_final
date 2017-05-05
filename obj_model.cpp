#include "obj_model.h"

using namespace std;

ObjModel::ObjModel(GLfloat* i_points, GLint* i_faces, GLfloat* i_tex_coords,  
                   int point_count, int face_count, int tex_count) {
  // Set defaults
  rot_speed = vec3(0.0, 0.0, 0.0);
  rot = vec3(1.0, 1.0, 1.0);
  obj_scale = vec3(1.0, 1.0, 1.0);
  translation = vec3(0, 0, 0);

  // Convert from arrays to vectors
  points = convert_to_vector<GLfloat>(i_points, point_count);
  faces = convert_to_vector<GLint>(i_faces, face_count);
  texture_coords = convert_to_vector<GLfloat>(i_tex_coords, tex_count);

  shader = 1;
}

/**
 * Helper function to convert from an array to a vector
 */
template <typename T> vector<T> ObjModel::convert_to_vector(T* array, int count) {
  // Create an empty vector
  vector<T> new_vec = { };

  for (int i = 0; i < count; i++) {
    new_vec.push_back(array[i]);
  }

  return new_vec;
}

template <typename T> T* ObjModel::convert_to_array(vector<T> vector) {
  T* array = new T[vector.size()];

  for (unsigned int i = 0; i < vector.size(); i++) {
    array[i] = vector.at(i);
  }
  
  return array;
}

void ObjModel::update(double seconds) {
  model_mat = rotate_x_deg(identity_mat4(), seconds * rot_speed.v[0]) *
              rotate_y_deg(identity_mat4(), seconds * rot_speed.v[1]) *
              rotate_z_deg(identity_mat4(), seconds * rot_speed.v[2]) *
              translate(identity_mat4(), translation) *
              scale(identity_mat4(), obj_scale);
}

// Getter methods
GLfloat* ObjModel::get_points() {
  return convert_to_array<GLfloat>(points);
}

GLint* ObjModel::get_faces() {
  return convert_to_array<GLint>(faces);
}

GLfloat* ObjModel::get_texture_coords() {
  return convert_to_array<GLfloat>(texture_coords);
}

mat4 ObjModel::get_model_mat() {
  return model_mat;
}

GLuint ObjModel::get_shader() {
  return shader;
}

int ObjModel::get_mat_location() {
  return mat_location;
}

void ObjModel::set_rot(float r_x, float r_y, float r_z) {
  rot = vec3(r_x, r_y, r_z);
}

void ObjModel::set_rot_speed(float speed_x, float speed_y, float speed_z) {
  rot_speed = vec3(speed_x, speed_y, speed_z);
}

void ObjModel::set_scale(double s_x, double s_y, double s_z) {
  obj_scale = vec3(s_x, s_y, s_z);
}

void ObjModel::set_translation(float x, float y, float z) {
  translation = vec3(x, y, z);
}

void ObjModel::set_shader(GLuint shader) {
  this->shader = shader;
}

void ObjModel::set_mat_location(int mat_location) {
  this->mat_location = mat_location;
}
