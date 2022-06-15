// This example is heavily based on the tutorial at https://open.gl
#include <vector>
#include <cmath>
#include <string>
// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>

// Linear Algebra Library
#include <Eigen/Core>

// Timer
#include <chrono>
using namespace Eigen;
using namespace std;

Program program;

struct Triangle {
	vector<Vector3d> Color;
	vector<Vector3d> vertices;
};

enum Mode {
	NONE,
	INSERTION,
	TRANSLATION,
	DELETION,
	COLOR_CHANGE,
	ROTATIONCLOCKWISE,
	ROTATIONCOUNTERCLOCK,
	SCALEDOWN,
	SCALEUP,
	ZOOMIN,
	ZOOMOUT,
	VIEWCONTROLW,
	VIEWCONTROLA,
	VIEWCONTROLS,
	VIEWCONTROLD
};

const double pi = 3.1415926;



// VertexBufferObject wrapper
VertexBufferObject triangleVBO;
VertexBufferObject lineVBO;
VertexBufferObject pointVBO;
VertexBufferObject triangleColorVBO;
VertexBufferObject lineColorVBO;
VertexBufferObject pointColorVBO;
VertexArrayObject VAO;
MatrixXf Data;
Mode mode = INSERTION;
typedef vector<Triangle> TriangleListType;
TriangleListType TriangleList;
vector <vector<Triangle>> Keyframes;
double a = .50;
int KeyframeIndex = 0;

int translating_triangle = -1;
// We also store the start of the translation
Vector2d translation_start;
float Zoom = 1;
Vector2d pan(0, 0);

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



double lerp(double x, double y, double a)
{
	return (x*(1 - a) + y * a);
}
Vector2d lerp(Vector2d v1, Vector2d v2, double a)
{
	Vector2d v;
	v.x() = lerp(v1.x(), v2.x(), a);
	v.y() = lerp(v1.y(), v2.y(), a);
	return v;
}

Vector3d lerp(Vector3d v1, Vector3d v2, double a)
{
	Vector3d v;
	v.x() = lerp(v1.x(), v2.x(), a);
	v.y() = lerp(v1.y(), v2.y(), a);
	v.z() = lerp(v1.z(), v2.z(), a);
	return v;
}
Triangle lerp(Triangle  t1, Triangle  t2, double a)
{
	Triangle t;
	t.vertices.resize(3);
	t.Color.resize(3);
	t.vertices[0] = lerp(t1.vertices[0], t2.vertices[0], a);
	t.Color[0] = lerp(t1.Color[0], t2.Color[0], a);
	t.vertices[1] = lerp(t1.vertices[1], t2.vertices[1], a);
	t.Color[1] = lerp(t1.Color[1], t2.Color[1], a);
	t.vertices[2] = lerp(t1.vertices[2], t2.vertices[2], a);
	t.Color[2] = lerp(t1.Color[2], t2.Color[2], a);
	return t;
}
Triangle MultipleKeyframeSupport(vector<Triangle> Keyframes,double a, int index =0)
{
	Triangle InterpolateValue;
	double Interpolation = a * (Keyframes.size()-1);
	InterpolateValue = lerp(Keyframes[index], Keyframes[index++], Interpolation);
	return InterpolateValue;
}
TriangleListType lerp(const TriangleListType& t1, const TriangleListType& t2, double a)
{
	TriangleListType output;
	// If the lists don't have the same size interpolate between the n first triangles and handle the remaining separately
	int size = std::min(t1.size(), t2.size());
	for (int i = 0; i < size; ++i)
	{
		output.push_back(lerp(t1[i], t2[i], a));
	}
	// handle the case where the size is different
	if (a == 0) {
		for (int i = size; i < t1.size(); ++i)
		{
			output.push_back(t1[i]);
		}
	}
	else if (a == 1) {
		for (int i = size; i < t2.size(); ++i)
		{
			output.push_back(t2[i]);
		}
	}
	return output;
}




Vector2d getMousePosition(GLFWwindow* window) {
	
	int width, height;
	glfwGetWindowSize(window, &width, &height);

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

// Given 2 vectors, this function computes the trnaslation vector from the start point to the end point
Vector2d translation(Vector2d start, Vector2d end) {
	
	return end - start;
}

// double area(int x1, int y1, int x2, int y2, int x3, int y3)
double area(double x1, double y1, double x2, double y2, double x3, double y3)
{
	return abs((x1*(y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
}


// Given a triangle and a point representing the mouse position, return true if the point is inside the triangle, false otherwise
bool isInsideTriangle(Vector2d point, const Triangle& triangle)
{
	if (triangle.vertices.size() < 3)
	{
		
		return false;
	}
	auto X1 = triangle.vertices[0].x();
	auto X2 = triangle.vertices[1].x();
	auto X3 = triangle.vertices[2].x();
	auto Y1 = triangle.vertices[0].y();
	auto Y2 = triangle.vertices[1].y();
	auto Y3 = triangle.vertices[2].y();
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
// Given a list of triangles and a point, return the index of the first triangle containing the point. Return -1 if no triangle contains the point
int isInsideTriangle(Vector2d point, const vector<Triangle>& triangleList) {

	for (int i = 0; i < triangleList.size(); ++i)
	{
		if (isInsideTriangle(point, triangleList[i]))
		{
			return i;
		}
	}
	return -1;
	
}

void scaleAndRotate(Triangle& triangle, double angle = 0, double scale = 1) {
	const double cosine = std::cos(angle);
	const double sine = std::sin(angle);

	Vector3d center = (triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2]) / 3.0;

	double X1 = triangle.vertices[0].x() - center.x();
	double X2 = triangle.vertices[1].x() - center.x();
	double X3 = triangle.vertices[2].x() - center.x();
	double Y1 = triangle.vertices[0].y() - center.y();
	double Y2 = triangle.vertices[1].y() - center.y();
	double Y3 = triangle.vertices[2].y() - center.y();

	// Apply rotation
	auto NX1 = X1 * cosine - Y1 * sine;
	auto NY1 = X1 * sine + Y1 * cosine;
	// Apply scale
	NX1 *= scale;
	NY1 *= scale;
	// Store result
	triangle.vertices[0].x() = NX1 + center.x();
	triangle.vertices[0].y() = NY1 + center.y();

	auto NX2 = X2 * cosine - Y2 * sine;
	auto NY2 = X2 * sine + Y2 * cosine;
	// Apply scale
	NX2 *= scale;
	NY2 *= scale;
	triangle.vertices[1].x() = NX2 + center.x();
	triangle.vertices[1].y() = NY2 + center.y();
	auto NX3 = X3 * cosine - Y3 * sine;
	auto NY3 = X3 * sine + Y3 * cosine;
	// Apply scale
	NX3 *= scale;
	NY3 *= scale;
	triangle.vertices[2].x() = NX3 + center.x();
	triangle.vertices[2].y() = NY3 + center.y();
}


void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	Vector2d mouse_pos = getMousePosition(window);
	// Update the position of the first vertex if the left button is pressed
	if (button == GLFW_MOUSE_BUTTON_LEFT)
	{
		if (action == GLFW_PRESS)
		{
			if (mode == INSERTION) {
				if (!TriangleList.empty() && TriangleList.back().vertices.size() < 3)
				{
					TriangleList.back().Color.push_back(Vector3d(1, 0, 0));
					TriangleList.back().vertices.push_back(Vector3d(mouse_pos.x(), mouse_pos.y(), 0));

				}
				else
				{
					Triangle New_Triangle;
					New_Triangle.Color.push_back(Vector3d(1, 0, 0));
					New_Triangle.vertices.push_back(Vector3d(mouse_pos.x(), mouse_pos.y(), 0));
					TriangleList.push_back(New_Triangle);
				}
			}
			else if (mode == TRANSLATION)
			{
				
				translating_triangle = isInsideTriangle(mouse_pos, TriangleList);
				translation_start = mouse_pos;
			}
			else if (mode == DELETION)
			{
				for (int i = 0; i < TriangleList.size(); ++i)
				{
					if (isInsideTriangle(mouse_pos, TriangleList[i]))
					{
						TriangleList.erase(TriangleList.begin() + i);
					}
				}
			}
			else if (mode == ROTATIONCLOCKWISE)
			{
				for (int i = 0; i < TriangleList.size(); ++i)
				{
					if (isInsideTriangle(mouse_pos, TriangleList[i]))
					{
						const double angle = (pi / 180) * 10;
						scaleAndRotate(TriangleList[i], angle);
					}
				}

			}
			else if (mode == SCALEUP)
			{
				for (int i = 0; i < TriangleList.size(); ++i)
				{
					if (isInsideTriangle(mouse_pos, TriangleList[i]))
					{
						const double angle = 0;
						const double scale_factor = 1.25;
						scaleAndRotate(TriangleList[i], angle, scale_factor);
					}
				}
			}
			else if (mode == SCALEDOWN)
			{
				for (int i = 0; i < TriangleList.size(); ++i)
				{
					if (isInsideTriangle(mouse_pos, TriangleList[i]))
					{
						const double angle = 0;
						const double scale_factor = 0.75;
						scaleAndRotate(TriangleList[i], angle, scale_factor);
					}
				}
			}
			else if (mode == ROTATIONCOUNTERCLOCK)
			{
				for (int i = 0; i < TriangleList.size(); ++i)
				{
					if (isInsideTriangle(mouse_pos, TriangleList[i]))
					{
						const double angle = -(pi / 180) * 10;
						scaleAndRotate(TriangleList[i], angle);
					}
				}
			}
			else if (mode == COLOR_CHANGE)
			{
				int triangle_ixd = -1;
				int vertex_ixd = -1;
				double closest_distance_squared = 9999999999.0;
				for (int i = 0; i < TriangleList.size(); ++i) {
					for (int j = 0; j < TriangleList[i].vertices.size(); ++j) {
						double vx = TriangleList[i].vertices[j].x();
						double vy = TriangleList[i].vertices[j].y();
						double mx = mouse_pos.x();
						double my = mouse_pos.y();
						double distance_squared = (vx - mx) * (vx - mx) + (vy - my) * (vy - my);
						if (distance_squared < closest_distance_squared) {
							closest_distance_squared = distance_squared;
							triangle_ixd = i;
							vertex_ixd = j;
						}
					}
				}
				// if the mouse distance to the vertex is more than a threshold (far away) then we don't pick the vertex
				const double distance_threshold = 0.02;
				if (closest_distance_squared < distance_threshold) {
					TriangleList[triangle_ixd].Color[vertex_ixd] = colorTable[currentColorIdx];
				}

			}
		}
		if (action == GLFW_RELEASE) {
			if (mode == TRANSLATION)
			{
				
				if (translating_triangle >= 0) {
					Vector2d t = translation(translation_start, mouse_pos);
					for (int i = 0; i < 3; ++i) {
						TriangleList[translating_triangle].vertices[i].x() += t.x();
						TriangleList[translating_triangle].vertices[i].y() += t.y();
					}
					translating_triangle = -1; // done translating
				}
			}
		}
	}
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	switch (key)
	{
	case  GLFW_KEY_I:
		mode = INSERTION;
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
	case  GLFW_KEY_W:
		mode = VIEWCONTROLW;
		pan += Vector2d(0, 0.20);
		break;
	case  GLFW_KEY_A:
		mode = VIEWCONTROLA;
		pan += Vector2d(-0.20, 0);
		break;
	case  GLFW_KEY_S:
		mode = VIEWCONTROLS;
		pan += Vector2d(0, -0.20);
		break;
	case  GLFW_KEY_D:
		mode = VIEWCONTROLD;
		pan += Vector2d(0.20, 0);
		break;
	case GLFW_KEY_KP_ADD:
		mode = ZOOMIN;
		break;
	case GLFW_KEY_KP_SUBTRACT:
		mode = ZOOMOUT;
		break;
	case GLFW_KEY_MINUS:
		Zoom *= 0.8;
		break;
	case GLFW_KEY_EQUAL:
		Zoom *= 1.2;
		break;
	case GLFW_KEY_KP_1: case GLFW_KEY_1:
		currentColorIdx = 0;
		break;
	case GLFW_KEY_KP_2: case GLFW_KEY_2:
		currentColorIdx = 1;
		break;
	case GLFW_KEY_KP_3: case GLFW_KEY_3:
		currentColorIdx = 2;
		break;
	case GLFW_KEY_KP_4: case GLFW_KEY_4:
		currentColorIdx = 3;
		break;
	case GLFW_KEY_KP_5: case GLFW_KEY_5:
		currentColorIdx = 4;
		break;
	case GLFW_KEY_KP_6: case GLFW_KEY_6:
		currentColorIdx = 5;
		break;
	case GLFW_KEY_KP_7: case GLFW_KEY_7:
		currentColorIdx = 6;
		break;
	case GLFW_KEY_KP_8: case GLFW_KEY_8:
		currentColorIdx = 7;
		break;
	case GLFW_KEY_KP_9: case GLFW_KEY_9:
		currentColorIdx = 8;
		break;
	case GLFW_KEY_SPACE:
		if (action == GLFW_RELEASE) {
			Keyframes.push_back(TriangleList);
		}
		break;
	case GLFW_KEY_LEFT:
		if (Keyframes.size() > 1)
		{
			a -= 0.02;
			a = max(a, 0.0);
			TriangleList = lerp(Keyframes[KeyframeIndex], Keyframes[KeyframeIndex++], a);
			KeyframeIndex++;
		}
		break;
	case GLFW_KEY_RIGHT:
		if (Keyframes.size() > 1)
		{
			a += 0.02;
			a = min(a, 1.0);
			TriangleList = lerp(Keyframes[KeyframeIndex], Keyframes[KeyframeIndex++], a);
			KeyframeIndex++;
		}
	
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
	glewExperimental = GL_TRUE;
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


	// Initialize the VBO with the vertices data
	// A VBO is a data container that lives in the GPU memory
	triangleVBO.init();
	lineVBO.init();
	pointVBO.init();

	triangleColorVBO.init();
	lineColorVBO.init();
	pointColorVBO.init();

	// Initialize the OpenGL Program
	// A program controls the OpenGL pipeline and it must contains
	// at least a vertex shader and a fragment shader to be valid
	Program program;
	const GLchar* vertex_shader =
		"#version 150 core\n"
		"in vec2 position;"
		"in vec3 vs_in_color;"
		"out vec3 vs_out_color;"
		"uniform float zoom;" // scale
		"uniform vec2 pan;" // translation
		"void main()"
		"{"
		"    gl_Position = vec4(position * zoom + pan, 0.0, 1.0);"
		"    vs_out_color = vs_in_color;"
		"}";
	const GLchar* fragment_shader =
		"#version 150 core\n"
		"in vec3 vs_out_color;"
		"out vec4 outColor;"
		"void main()"
		"{"
		"    outColor = vec4(vs_out_color, 1.0);"
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

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window))
	{
		int vw, vh;
		glfwGetFramebufferSize(window, &vw, &vh);
		glViewport(0, 0, vw, vh);
		int triangle_count = 0;
		int point_count = 0;
		int lines_count = 0;
		for (int i = 0; i < TriangleList.size(); ++i) {
			if (TriangleList[i].vertices.size() == 2)
			{
				++lines_count;
			}
			else if (TriangleList[i].vertices.size() == 1) {
				++point_count;
			}
			else if (TriangleList[i].vertices.size() == 3) {
				++triangle_count;
			}
		}


		MatrixXf triangle_matrix(2, triangle_count * 3);
		triangle_matrix.resize(2, triangle_count * 3);
		MatrixXf triangle_color_matrix(3, triangle_count * 3);
		triangle_color_matrix.resize(3, triangle_count * 3);

		MatrixXf line_matrix(2, lines_count * 2);
		line_matrix.resize(2, lines_count * 2);
		MatrixXf line_color_matrix(3, lines_count * 2);
		line_color_matrix.resize(3, lines_count * 2);

		MatrixXf point_matrix(2, point_count * 1);
		point_matrix.resize(2, point_count * 1);
		MatrixXf point_color_matrix(3, point_count * 1);
		point_color_matrix.resize(3, point_count * 1);

		int point_index = 0;
		int triangle_index = 0;
		int line_index = 0;
		for (int i = 0; i < TriangleList.size(); ++i) {
			if (TriangleList[i].vertices.size() == 2)
			{
				for (int j = 0; j < 2; ++j)
				{
					line_matrix(0, line_index * 2 + j) = TriangleList[i].vertices[j].x();
					line_matrix(1, line_index * 2 + j) = TriangleList[i].vertices[j].y();

					line_color_matrix(0, line_index * 2 + j) = 1;
					line_color_matrix(1, line_index * 2 + j) = 1;
					line_color_matrix(2, line_index * 2 + j) = 1;
				}
				++line_index;
			}
			else if (TriangleList[i].vertices.size() == 1) {
				for (int j = 0; j < 1; ++j)
				{
					point_matrix(0, point_index + j) = TriangleList[i].vertices[j].x();
					point_matrix(1, point_index + j) = TriangleList[i].vertices[j].y();

					point_color_matrix(0, point_index + j) = 1;
					point_color_matrix(1, point_index + j) = 1;
					point_color_matrix(2, point_index + j) = 1;
				}
				++point_index;
			}
			else if (TriangleList[i].vertices.size() == 3) {
				Vector2d t = Vector2d(0, 0);
				if (i == translating_triangle) {
					
					// The triangle is being translating
					// Apply the translation
					Vector2d mouse_pos = getMousePosition(window);
					t = translation(translation_start, mouse_pos);
				}
				for (int j = 0; j < 3; ++j)
				{
					triangle_matrix(0, triangle_index * 3 + j) = TriangleList[i].vertices[j].x() + t.x();
					triangle_matrix(1, triangle_index * 3 + j) = TriangleList[i].vertices[j].y() + t.y();

					triangle_color_matrix(0, triangle_index * 3 + j) = TriangleList[i].Color[j].x();
					triangle_color_matrix(1, triangle_index * 3 + j) = TriangleList[i].Color[j].y();
					triangle_color_matrix(2, triangle_index * 3 + j) = TriangleList[i].Color[j].z();
				}
				++triangle_index;
			}
		}
		triangleVBO.update(triangle_matrix);
		lineVBO.update(line_matrix);
		pointVBO.update(point_matrix);

		triangleColorVBO.update(triangle_color_matrix);
		lineColorVBO.update(line_color_matrix);
		pointColorVBO.update(point_color_matrix);

		// Bind your VAO (not necessary if you have only one)
		VAO.bind();

		// Bind your program
		program.bind();


		// Clear the framebuffer
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		glUniform1f(program.uniform("zoom"), Zoom); // a scale of 1 means no zoom
		// TODO pass the translation to the vertex shader
		glUniform2f(program.uniform("pan"), pan.x(), pan.y()); // a zero translation means no pan

		// glUniform3f(program.uniform("triangleColor"), 1, 0, 0);
		if (triangle_count > 0)
		{
			program.bindVertexAttribArray("position", triangleVBO);
			program.bindVertexAttribArray("vs_in_color", triangleColorVBO);
			glDrawArrays(GL_TRIANGLES, 0, triangle_count * 3);
		}

		// glUniform3f(program.uniform("triangleColor"), 1, 1, 1);
		if (lines_count > 0)
		{
			program.bindVertexAttribArray("position", lineVBO);
			program.bindVertexAttribArray("vs_in_color", lineColorVBO);
			glDrawArrays(GL_LINES, 0, lines_count * 2);
		}

		if (point_count > 0)
		{
			program.bindVertexAttribArray("position", pointVBO);
			program.bindVertexAttribArray("vs_in_color", pointColorVBO);
			glDrawArrays(GL_POINTS, 0, point_count);
		}

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();
	}

	// Deallocate opengl memory
	program.free();
	VAO.free();
	triangleVBO.free();

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}
