
// This example is heavily based on the tutorial at https://open.gl

#include <vector>
#include <cmath>
#include <string>
// OpenGL Helpers to reduce the clutter
#include "Helpers.h"
#include <fstream>
#include <iostream>
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>

// Linear Algebra Library
#include <Eigen/Core>
#include <Eigen/Geometry>

// Timer
#include <chrono>
using namespace Eigen;
using namespace std;
MatrixXf projection(4, 4);
MatrixXf camera(4, 4);
struct Triangle {
	Vector3d vertices[3];
	Vector3d normal; // geometric triangle normal
	Vector3d smooth_normals[3]; //  per-vertex averaged normal
};

VertexArrayObject VAO;

typedef vector<Triangle> TriangleListType;

struct Object
{
	TriangleListType triangles;
	Transform<float, 3, Affine> transform;
	VertexBufferObject VBO; // vertex buffer object (vertices only)
	VertexBufferObject NBO; // normal buffer object (normal)
	Vector3d color;
};

// Create a list of object
std::vector<Object> Objects;

// Contains the vertex positions
Vector3d colorTable[9] = {
	Vector3d(1,0,0),
	Vector3d(0,1,0),
	Vector3d(0,0,1),
	Vector3d(0,1,1),
	Vector3d(1,1,0),
	Vector3d(1,0,1),
	Vector3d(0.5,0.25,0.75),
	Vector3d(0.5,0.75,0.25),
	Vector3d(0.75,0.5,0.25),
};
int currentColorIdx = 0;

float Zoom = 1;
Vector2d pan(0, 0);

// Same as the 2d editor
int translating_triangle = -1;
Vector2d translation_start;

int rotation_axis = 0; // 0 means rotate around x, 1 around y, 2 around z

float camera_angle = 0;
float camera_elevation = 0;
float camera_distance = -3;

enum Mode {
	NONE,
	UNITCUBE,
	BUMPYCUBE,
	BUNNY,
	TRANSLATION,
	DELETION,
	COLOR_CHANGE,
	ROTATIONCLOCKWISE,
	ROTATIONCOUNTERCLOCK,
	SCALEDOWN,
	SCALEUP,
	WIREFRAME,
	FLATSHADING,
	PHONGSHADING,
	ORTHOGRAPHIC,
	PERSPECTIVE
};

Mode mode = UNITCUBE;
Mode lighting_mode = WIREFRAME;
Mode projection_mode = ORTHOGRAPHIC;

const float canvas_size = 3;


void setProjection(float fovY, float aspectRatio, float N, float F, float r = -canvas_size, float l = canvas_size, float t = canvas_size, float b = -canvas_size)
{
	projection = MatrixXf::Identity(4, 4);
	if (projection_mode == PERSPECTIVE) {

		projection <<
			2 * N / (r - l), 0, (r + l) / (r - l), 0,
			0, 2 * N / (t - b), (t + b) / (t - b), 0,
			0, 0, (F + N) / (N - F), (2 * N*F) / (N - F),
			0, 0, -1, 0;
	}
	else {
		projection <<
			2 / (r - l), 0, 0, 0,
			0, 2 / (t - b), 0, 0,
			0, 0, -2 / (F - N), 0,
			-(r + l) / (r - l), -(t + b) / (t - b), -(F + N) / (F - N), 1;
	}
}

void makeCameraMatrix() {
	Vector3d camera_position;
	camera_position.z() = cos(camera_angle) * camera_distance;
	camera_position.x() = sin(camera_angle) * camera_distance;
	camera_position.y() = camera_elevation;

	Vector3d world_origin = Vector3d(0, 0, 0);

	Vector3d up_vector = Vector3d(0, 1, 0);
	Vector3d forward_vector = -(world_origin - camera_position).normalized();
	if (world_origin == camera_position) {
		forward_vector = Vector3d(0, 0, 1);
	}
	Vector3d right_vector = up_vector.cross(forward_vector);
	up_vector = forward_vector.cross(right_vector);

	/*
	camera << right_vector.x(), right_vector.y(), right_vector.z(), -right_vector.dot(camera_position),
	  up_vector.x(), up_vector.y(), up_vector.z(), -up_vector.dot(camera_position),
	  forward_vector.x(), forward_vector.y(), forward_vector.z(), -forward_vector.dot(camera_position),
	  0, 0, 0, 1;
	*/
	camera << right_vector.x(), up_vector.x(), forward_vector.x(), -right_vector.dot(camera_position),
		right_vector.y(), up_vector.y(), forward_vector.y(), -up_vector.dot(camera_position),
		right_vector.z(), up_vector.z(), forward_vector.z(), -forward_vector.dot(camera_position),
		0, 0, 0, 1;
}

void generateNormal(Triangle& triangle) {
	//ABC triangle
	Vector3d A = triangle.vertices[0];
	Vector3d B = triangle.vertices[1];
	Vector3d C = triangle.vertices[2];
	Vector3d V = B - A;
	Vector3d W = C - A;
	Vector3d Normal = V.cross(W);
	if (Normal.x() != 0 || Normal.y() != 0 || Normal.z() != 0) {
		triangle.normal = Normal.normalized();
	}
	else {
		triangle.normal = Normal;
	}
}

void generateNormals(TriangleListType& triangles) {
	// fill the normal
	for (int i = 0; i < triangles.size(); ++i)
	{
		generateNormal(triangles[i]);
	}
}

TriangleListType loadOff(const char* path) {
	TriangleListType triangles;

	std::ifstream input_file(path);
	if (!input_file.is_open()) {
		cout << "ERROR can't open " << path << "\n";
		return triangles;
	}

	string identifier;
	input_file >> identifier;
	if (identifier != "OFF") {
		cout << "ERROR missing OFF\n";
		return triangles;
	}
	int vert_count = 0;
	int face_count = 0;
	int other = 0; // TODO
	input_file >> vert_count >> face_count >> other;

	std::vector<Vector3d> verts;
	Vector3d average = Vector3d(0, 0, 0);
	for (int i = 0; i < vert_count; ++i) {
		float v[3];
		input_file >> v[0] >> v[1] >> v[2];
		verts.push_back(Vector3d(v[0], v[1], v[2]));
		average += verts.back();
	}
	// Center it
	average = average / double(vert_count);
	for (int i = 0; i < vert_count; ++i) {
		verts[i] = verts[i] - average;
	}
	// Scale it to unit cube
	double max_dist = 0;
	for (int i = 0; i < vert_count; ++i) {
		double d = verts[i].x() * verts[i].x() + verts[i].y() * verts[i].y() + verts[i].z() * verts[i].z();
		max_dist = max(d, max_dist);
	}
	if (max_dist > 0) {
		max_dist = sqrt(max_dist);
		for (int i = 0; i < vert_count; ++i) {
			verts[i] = verts[i] / max_dist;
		}
	}

	triangles.resize(face_count);

	// Initialize smooth_normals to 0
	std::vector<Vector3d> smooth_normals(verts.size(), Vector3d(0, 0, 0));

	// Read faces certex indices
	std::vector<int> indices;
	for (int i = 0; i < face_count; ++i) {
		int face_vert_count;
		input_file >> face_vert_count;
		if (face_vert_count != 3) {
			cout << "ERROR faces must be triangles\n";
			return triangles;
		}
		int v1, v2, v3;
		input_file >> v1 >> v2 >> v3;
		indices.push_back(v1);
		indices.push_back(v2);
		indices.push_back(v3);
	}

	for (int i = 0; i < face_count; ++i) {
		int v1 = indices[i * 3];
		int v2 = indices[i * 3 + 1];
		int v3 = indices[i * 3 + 2];
		triangles[i].vertices[0] = verts[v1];
		triangles[i].vertices[1] = verts[v2];
		triangles[i].vertices[2] = verts[v3];
		generateNormal(triangles[i]);

		// Average normals for every vertex by summing them
		smooth_normals[v1] += triangles[i].normal;
		smooth_normals[v2] += triangles[i].normal;
		smooth_normals[v3] += triangles[i].normal;
	}

	// Save smooth normals and make sure they are normalized
	for (int i = 0; i < face_count; ++i) {
		int v1 = indices[i * 3];
		int v2 = indices[i * 3 + 1];
		int v3 = indices[i * 3 + 2];
		triangles[i].smooth_normals[0] = smooth_normals[v1].normalized();
		triangles[i].smooth_normals[1] = smooth_normals[v2].normalized();
		triangles[i].smooth_normals[2] = smooth_normals[v3].normalized();
	}

	return triangles;
}

TriangleListType makeCube() {
	TriangleListType CubeTriangleList;
	Triangle t;
	// We have 6 faces

	// Front face
	t.vertices[0].x() = -1; t.vertices[0].y() = -1; t.vertices[0].z() = 1;
	t.vertices[1].x() = 1; t.vertices[1].y() = -1; t.vertices[1].z() = 1;
	t.vertices[2].x() = 1; t.vertices[2].y() = 1; t.vertices[2].z() = 1;
	CubeTriangleList.push_back(t);
	t.vertices[0].x() = 1; t.vertices[0].y() = 1; t.vertices[0].z() = 1;
	t.vertices[1].x() = -1; t.vertices[1].y() = 1; t.vertices[1].z() = 1;
	t.vertices[2].x() = -1; t.vertices[2].y() = -1; t.vertices[2].z() = 1;
	CubeTriangleList.push_back(t);

	// Back face
	t.vertices[0].x() = 1; t.vertices[0].y() = -1; t.vertices[0].z() = -1;
	t.vertices[1].x() = -1; t.vertices[1].y() = -1; t.vertices[1].z() = -1;
	t.vertices[2].x() = -1; t.vertices[2].y() = 1; t.vertices[2].z() = -1;
	CubeTriangleList.push_back(t);
	t.vertices[0].x() = -1; t.vertices[0].y() = 1; t.vertices[0].z() = -1;
	t.vertices[1].x() = 1;  t.vertices[1].y() = 1; t.vertices[1].z() = -1;
	t.vertices[2].x() = 1; t.vertices[2].y() = -1; t.vertices[2].z() = -1;
	CubeTriangleList.push_back(t);

	// Top face
	t.vertices[0].x() = -1; t.vertices[0].y() = 1; t.vertices[0].z() = 1;
	t.vertices[1].x() = 1; t.vertices[1].y() = 1; t.vertices[1].z() = 1;
	t.vertices[2].x() = 1; t.vertices[2].y() = 1; t.vertices[2].z() = -1;
	CubeTriangleList.push_back(t);
	t.vertices[0].x() = 1; t.vertices[0].y() = 1; t.vertices[0].z() = -1;
	t.vertices[1].x() = -1; t.vertices[1].y() = 1; t.vertices[1].z() = -1;
	t.vertices[2].x() = -1; t.vertices[2].y() = 1; t.vertices[2].z() = 1;
	CubeTriangleList.push_back(t);

	// Bottom face
	t.vertices[0].x() = -1; t.vertices[0].y() = -1; t.vertices[0].z() = -1;
	t.vertices[1].x() = 1; t.vertices[1].y() = -1; t.vertices[1].z() = -1;
	t.vertices[2].x() = 1; t.vertices[2].y() = -1; t.vertices[2].z() = 1;
	CubeTriangleList.push_back(t);
	t.vertices[0].x() = 1; t.vertices[0].y() = -1; t.vertices[0].z() = 1;
	t.vertices[1].x() = -1; t.vertices[1].y() = -1; t.vertices[1].z() = 1;
	t.vertices[2].x() = -1; t.vertices[2].y() = -1; t.vertices[2].z() = -1;
	CubeTriangleList.push_back(t);

	// left side
	t.vertices[0].x() = -1; t.vertices[0].y() = -1; t.vertices[0].z() = -1;
	t.vertices[1].x() = -1; t.vertices[1].y() = -1; t.vertices[1].z() = 1;
	t.vertices[2].x() = -1; t.vertices[2].y() = 1; t.vertices[2].z() = 1;
	CubeTriangleList.push_back(t);
	t.vertices[0].x() = -1; t.vertices[0].y() = 1; t.vertices[0].z() = 1;
	t.vertices[1].x() = -1; t.vertices[1].y() = 1; t.vertices[1].z() = -1;
	t.vertices[2].x() = -1; t.vertices[2].y() = -1; t.vertices[2].z() = -1;
	CubeTriangleList.push_back(t);

	// right side
	t.vertices[0].x() = 1; t.vertices[0].y() = -1; t.vertices[0].z() = 1;
	t.vertices[1].x() = 1; t.vertices[1].y() = -1; t.vertices[1].z() = -1;
	t.vertices[2].x() = 1; t.vertices[2].y() = 1; t.vertices[2].z() = -1;
	CubeTriangleList.push_back(t);
	t.vertices[0].x() = 1; t.vertices[0].y() = 1; t.vertices[0].z() = -1;
	t.vertices[1].x() = 1; t.vertices[1].y() = 1; t.vertices[1].z() = 1;
	t.vertices[2].x() = 1; t.vertices[2].y() = -1; t.vertices[2].z() = 1;
	CubeTriangleList.push_back(t);

	generateNormals(CubeTriangleList);

	
   // Add smooth normals
	for (int i = 0; i < CubeTriangleList.size(); ++i) {
		CubeTriangleList[i].smooth_normals[0] = CubeTriangleList[i].normal;
		CubeTriangleList[i].smooth_normals[1] = CubeTriangleList[i].normal;
		CubeTriangleList[i].smooth_normals[2] = CubeTriangleList[i].normal;
	}

	return CubeTriangleList;
}

// double area(int x1, int y1, int x2, int y2, int x3, int y3)
double area(double x1, double y1, double x2, double y2, double x3, double y3)
{
	return abs((x1*(y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
}

// Given a triangle and a point representing the mouse position, return true if the point is inside the triangle, false otherwise
bool isInsideTriangle(Vector2d point, const Triangle& triangle, Transform<float, 3, Affine> transform)
{
	Vector2d transformed[3];
	for (int i = 0; i < 3; ++i) {
		MatrixXf model_view_proj = projection * camera * transform.matrix();
		Vector4f v = model_view_proj * Vector4f(triangle.vertices[i].x(), triangle.vertices[i].y(), triangle.vertices[i].z(), 1.f);
		transformed[i] = Vector2d(v.x() / v.w(), v.y() / v.w());
	}

	auto X1 = transformed[0].x();
	auto X2 = transformed[1].x();
	auto X3 = transformed[2].x();
	auto Y1 = transformed[0].y();
	auto Y2 = transformed[1].y();
	auto Y3 = transformed[2].y();
	auto A = area(X1, Y1, X2, Y2, X3, Y3);
	auto XPoint = point.x();
	auto YPoint = point.y();
	float A1 = area(XPoint, YPoint, X2, Y2, X3, Y3);

	float A2 = area(X1, Y1, XPoint, YPoint, X3, Y3);

	float A3 = area(X1, Y1, X2, Y2, XPoint, YPoint);

	/* Check if sum of A1, A2 and A3 is same as A */


	const double epsilon = 0.0001;
	return abs(A - (A1 + A2 + A3)) < epsilon;

}

bool isInsideTriangle(Vector2d point, const Object& object) {
	for (int i = 0; i < object.triangles.size(); ++i) {
		if (isInsideTriangle(point, object.triangles[i], object.transform)) {
			return true;
		}
	}
	return false;
}

// return the index of the object containing the triangle
int isInsideTriangle(Vector2d point) {
	for (int i = 0; i < Objects.size(); ++i) {
		if (isInsideTriangle(point, Objects[i])) {
			return i;
		}
	}
	return -1;
}

Vector2d getMousePosition(GLFWwindow* window) {
	


	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	float aspectRatio = float(height) / float(width);

	// Get the position of the mouse in the window
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	// Convert screen position to world coordinates
	double xworld = ((xpos / double(width)) * 2) - 1;
	double yworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw

	
	// new_position = position * zoom + pan
	// Let's invert the function (step by step)
	// new_position = position * zoom + pan
	// new_position - pan = position * zoom 
	// (new_position - pan) / zoom = position
	// So I you can see, you have to first substract pan and then divide by zoom
	// Before you were cancelling the zoom and the pan

	// Cancel the translation byy applying the inverse (not a division. It's substraction)
	xworld -= pan.x();
	yworld -= pan.y();

	// Cancel the zoom by applying the inverse
	xworld /= Zoom;
	yworld /= Zoom;

	return Vector2d(xworld, yworld);
}

void updateNBO(Object& O) {
	MatrixXf normal_matrix(3, O.triangles.size() * 3);
	normal_matrix.resize(3, O.triangles.size() * 3);

	for (int i = 0; i < O.triangles.size(); ++i) {
		//generateNormlas(O.triangles);
		for (int j = 0; j < 3; ++j)
		{
			if (lighting_mode != PHONGSHADING) {
				normal_matrix(0, i * 3 + j) = O.triangles[i].normal.x();
				normal_matrix(1, i * 3 + j) = O.triangles[i].normal.y();
				normal_matrix(2, i * 3 + j) = O.triangles[i].normal.z();
			}
			else {
				normal_matrix(0, i * 3 + j) = O.triangles[i].smooth_normals[j].x();
				normal_matrix(1, i * 3 + j) = O.triangles[i].smooth_normals[j].y();
				normal_matrix(2, i * 3 + j) = O.triangles[i].smooth_normals[j].z();
			}
		}
	}
	O.NBO.update(normal_matrix);
}

void updateAllNBOs() {
	for (int i = 0; i < Objects.size(); ++i)
	{
		updateNBO(Objects[i]);
	}
}

void FillVBOObject(Object& O)
{
	MatrixXf triangle_matrix(3, O.triangles.size() * 3);
	triangle_matrix.resize(3, O.triangles.size() * 3);
	for (int i = 0; i < O.triangles.size(); ++i) {
		for (int j = 0; j < 3; ++j)
		{
			triangle_matrix(0, i * 3 + j) = O.triangles[i].vertices[j].x();
			triangle_matrix(1, i * 3 + j) = O.triangles[i].vertices[j].y();
			triangle_matrix(2, i * 3 + j) = O.triangles[i].vertices[j].z();
		}
	}
	O.VBO.update(triangle_matrix);

	updateNBO(O);
}


void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the position of the mouse in the window
	Vector2d mouse_pos = getMousePosition(window);
	// Update the position of the first vertex if the left button is pressed
	if (button == GLFW_MOUSE_BUTTON_LEFT)
	{
		const float rotation_angle = 0.1;
		if (action == GLFW_PRESS)
		{
			int idx = isInsideTriangle(mouse_pos);
			if (mode == TRANSLATION)
			{
				translating_triangle = idx;
				translation_start = mouse_pos;
			}
			else if (mode == COLOR_CHANGE)
			{
				if (idx >= 0) {
					Objects[idx].color = colorTable[currentColorIdx];
					currentColorIdx++;
					currentColorIdx = currentColorIdx % 9;
				}
			}
			else if (mode == SCALEUP)
			{
				if (idx >= 0) {
					Objects[idx].transform.scale(1.01);
				}
			}
			else if (mode == SCALEDOWN)
			{
				if (idx >= 0) {
					Objects[idx].transform.scale(0.99);
				}
			}
			else if (mode == DELETION)
			{
				if (idx >= 0) {
					swap(Objects[idx], Objects.back());
					Objects.pop_back();
				}
			}
			else if (mode == ROTATIONCOUNTERCLOCK)
			{
				if (idx >= 0) {
					float angle = 0.01;
					if (rotation_axis == 0) {
						Objects[idx].transform.rotate(AngleAxisf(rotation_angle, Vector3f::UnitX()));
					}
					if (rotation_axis == 1) {
						Objects[idx].transform.rotate(AngleAxisf(rotation_angle, Vector3f::UnitY()));
					}
					else if (rotation_axis == 2) {
						Objects[idx].transform.rotate(AngleAxisf(rotation_angle, Vector3f::UnitZ()));
					}
				}
			}
			else if (mode == ROTATIONCLOCKWISE)
			{
				if (idx >= 0) {
					float angle = 0.01;
					if (rotation_axis == 0) {
						Objects[idx].transform.rotate(AngleAxisf(-rotation_angle, Vector3f::UnitX()));
					}
					if (rotation_axis == 1) {
						Objects[idx].transform.rotate(AngleAxisf(-rotation_angle, Vector3f::UnitY()));
					}
					else if (rotation_axis == 2) {
						Objects[idx].transform.rotate(AngleAxisf(-rotation_angle, Vector3f::UnitZ()));
					}
				}
			}
		}
		if (action == GLFW_RELEASE) {
			if (mode == TRANSLATION)
			{
				if (translating_triangle >= 0) {
					Vector2d t = (mouse_pos - translation_start) * canvas_size * 2;
					Objects[translating_triangle].transform.translate(Vector3f(t.x(), t.y(), 0));
					translating_triangle = -1; // done translating
				}
			}
		}
	}
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	if (action == GLFW_PRESS) {
		switch (key)
		{
		case GLFW_KEY_SPACE:
			rotation_axis = (rotation_axis + 1) % 3; // 0 means rotate around x, 1 around y, 2 around z
			break;
		case GLFW_KEY_KP_1: case GLFW_KEY_1:
		{
			Object cube;
			cube.triangles = makeCube();
			cube.transform = cube.transform.Identity();
			cube.color = Vector3d(1, 1, 1);
			cube.VBO.init();
			cube.NBO.init();
			FillVBOObject(cube);
			Objects.push_back(cube);
		}
		break;
		case GLFW_KEY_KP_3: case GLFW_KEY_3:
		{
			Object bunny;
			bunny.triangles = loadOff(DATA_FOLDER "/bunny.off");
			bunny.transform = bunny.transform.Identity();
			bunny.color = Vector3d(1, 1, 1);
			bunny.VBO.init();
			bunny.NBO.init();
			FillVBOObject(bunny);
			Objects.push_back(bunny);
		}
		break;
		case GLFW_KEY_KP_2: case GLFW_KEY_2:
		{
			Object bumpy_cube;
			bumpy_cube.triangles = loadOff(DATA_FOLDER "/bumpy_cube.off");
			bumpy_cube.transform = bumpy_cube.transform.Identity();
			bumpy_cube.color = Vector3d(1, 1, 1);
			bumpy_cube.VBO.init();
			bumpy_cube.NBO.init();
			FillVBOObject(bumpy_cube);
			Objects.push_back(bumpy_cube);
		}
		break;
		case  GLFW_KEY_ESCAPE:
			mode = NONE;
			break;
		case  GLFW_KEY_O:
			mode = TRANSLATION;
			break;
		case  GLFW_KEY_P:
			mode = DELETION;
			break;
		case  GLFW_KEY_H:
			mode = ROTATIONCLOCKWISE;
			break;
		case  GLFW_KEY_J:
			mode = ROTATIONCOUNTERCLOCK;
			break;
		case  GLFW_KEY_K:
			mode = SCALEUP;
			break;
		case  GLFW_KEY_L:
			mode = SCALEDOWN;
			break;
		case  GLFW_KEY_C:
			mode = COLOR_CHANGE;
			break;
		case GLFW_KEY_F: //orthographic
			projection_mode = ORTHOGRAPHIC;
			break;
		case GLFW_KEY_G: //persepctive
			projection_mode = PERSPECTIVE;
			break;
		case GLFW_KEY_W: //wireframe
			lighting_mode = WIREFRAME;
			updateAllNBOs();
			break;
		case GLFW_KEY_E: //flat shading
			lighting_mode = FLATSHADING;
			updateAllNBOs();
			break;
		case GLFW_KEY_R: //phong shading
			lighting_mode = PHONGSHADING;
			updateAllNBOs();
			break;
		default:
			break;
		}
	}

	const float camera_angle_step = 0.1;
	switch (key) {
	case GLFW_KEY_RIGHT:
		camera_angle += camera_angle_step;
		makeCameraMatrix();
		break;
	case GLFW_KEY_LEFT:
		camera_angle -= camera_angle_step;
		makeCameraMatrix();
		break;
	default:
		break;
	}
}

int main(void)
{
	GLFWwindow* window;

	// Initialize the library
	if (!glfwInit())
		return -1;

	// Activate supersampling
	glfwWindowHint(GLFW_SAMPLES, 8);

	// Ensure that we get at least a 3.2 context
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

	// On apple we have to load a core profile with forward compatibility
#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

	// Create a windowed mode window and its OpenGL context
	window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	// Make the window's context current
	glfwMakeContextCurrent(window);

#ifndef __APPLE__
	glewExperimental = true;
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
#endif

	int major, minor, rev;
	major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
	minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
	rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
	printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
	printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
	printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

	// Initialize the VAO
	// A Vertex Array Object (or VAO) is an object that describes how the vertex
	// attributes are stored in a Vertex Buffer Object (or VBO). This means that
	// the VAO is not the actual object storing the vertex data,
	// but the descriptor of the vertex data.
	VertexArrayObject VAO;
	VAO.init();
	VAO.bind();

	/*
	Object cube;
	cube.triangles = makeCube();
	cube.transform = cube.transform.Identity();
	cube.transform.scale(0.5);
	*/


	// Initialize the OpenGL Program
	// A program controls the OpenGL pipeline and it must contains
	// at least a vertex shader and a fragment shader to be valid
	Program program;
	const GLchar* vertex_shader =
		"#version 150 core\n"
		"in vec3 position;"
		"in vec3 normal;"
		"out vec3 vs_normal;"
		"uniform float zoom;" // scale
		"uniform vec2 pan;" // translation
		"uniform mat4 model_view_proj;" // combined model, view and projection matrices
		"uniform mat4 model_view;" // combined model, view and projection matrices
		"void main()"
		"{"
		"    vs_normal = mat3(model_view) * normal;"
		"    vs_normal = normalize(vs_normal);"
		"    gl_Position = model_view_proj * vec4(position, 1.0);"
		"    gl_Position = vec4(gl_Position.xy * zoom + pan, gl_Position.zw);"
		"}";

	const GLchar* fragment_shader =
		"#version 150 core\n"
		"uniform vec4 color;" // the alpha of zero means light is off and 1 means lighting is on
		"out vec4 outColor;"
		"in vec3 vs_normal;"
		"void main()"
		"{"
		"    vec3 n = normalize(vs_normal);"
		"    if (color.a == 0) "
		"      outColor = vec4(0.5);"
		"    else "
		"      outColor = vec4(color.rgb * abs(n.z), 1.0);"
		"}";

	// Compile the two shaders and upload the binary to the GPU
	// Note that we have to explicitly specify that the output "slot" called outColor
	// is the one that we want in the fragment buffer (and thus on screen)
	program.init(vertex_shader, fragment_shader, "outColor");
	program.bind();

	// Save the current time --- it will be used to dynamically change the triangle color
	auto t_start = std::chrono::high_resolution_clock::now();

	// Register the keyboard callback
	glfwSetKeyCallback(window, key_callback);

	// Register the mouse callback
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	makeCameraMatrix();

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window))
	{
		glEnable(GL_DEPTH_TEST);

		// Update viewport
		int vw, vh;
		glfwGetFramebufferSize(window, &vw, &vh);
		glViewport(0, 0, vw, vh);

		// Bind your VAO (not necessary if you have only one)
		VAO.bind();

		// Bind your program
		program.bind();

		// Set the uniform value depending on the time difference
		auto t_now = std::chrono::high_resolution_clock::now();
		float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();

		// Clear the framebuffer
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// TODO an options to switch between wireframe and shaded triangle
		if (lighting_mode == WIREFRAME) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else {
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		glUniform1f(program.uniform("zoom"), Zoom); // a scale of 1 means no zoom
		// TODO pass the translation to the vertex shader
		glUniform2f(program.uniform("pan"), pan.x(), pan.y()); // a zero translation means no pan

		setProjection(45, float(vw) / vw, 0.1, 100);
		MatrixXf view_proj = projection * camera;

		// Draw
		// TODO put this in a loop for all objects
		for (int i = 0; i < Objects.size(); ++i)
		{
			MatrixXf model_view_proj = view_proj * Objects[i].transform.matrix();
			MatrixXf model_view = camera * Objects[i].transform.matrix();
			glUniformMatrix4fv(program.uniform("model_view_proj"), 1, false, model_view_proj.data());
			glUniformMatrix4fv(program.uniform("model_view"), 1, false, model_view.data());
			program.bindVertexAttribArray("normal", Objects[i].NBO);
			program.bindVertexAttribArray("position", Objects[i].VBO);
			glUniform4f(program.uniform("color"), Objects[i].color.x(), Objects[i].color.y(), Objects[i].color.z(), 1);
			glDrawArrays(GL_TRIANGLES, 0, Objects[i].triangles.size() * 3);
		}

		if (lighting_mode == FLATSHADING) {
			glDisable(GL_DEPTH_TEST);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			for (int i = 0; i < Objects.size(); ++i)
			{
				MatrixXf model_view_proj = view_proj * Objects[i].transform.matrix();
				MatrixXf model_view = camera * Objects[i].transform.matrix();
				glUniformMatrix4fv(program.uniform("model_view_proj"), 1, false, model_view_proj.data());
				glUniformMatrix4fv(program.uniform("model_view"), 1, false, model_view.data());
				program.bindVertexAttribArray("normal", Objects[i].NBO);
				program.bindVertexAttribArray("position", Objects[i].VBO);
				glUniform4f(program.uniform("color"), Objects[i].color.x(), Objects[i].color.y(), Objects[i].color.z(), 0);
				glDrawArrays(GL_TRIANGLES, 0, Objects[i].triangles.size() * 3);
			}
		}

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();
	}

	// Deallocate opengl memory
	program.free();
	VAO.free();
	for (int i = 0; i < Objects.size(); ++i) {
		Objects[i].VBO.free();
		Objects[i].NBO.free();
	}

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}
