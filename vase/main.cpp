#include <GL/glut.h>
#include <bevgrafmath2017.h>
#include <math.h>

#include <fstream>	

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

GLfloat perspectiveNumber = 100;
GLfloat pNumberShange = 0.1;
GLsizei winWidth = 800, winHeight = 600;

mat4 w2v,scaled, projection[2], translation[2];

float vaseRotate = pi() / 6;
GLfloat tInc = 0.05;

GLint keyStates[256];

float incTetha = 0.1;
float incU = 0.1;
float changeIncVals = 0.01;
float alfa = 1;

struct face {
	std::vector<vec3> p;
	vec3 normal;
};

//vase controllpoints
vec2 controlPoints[6] = {
	{1,0},	//p0
	{3,1},	//p1
	{3,2},	//p2
	{1,2},	//p3
	{1,3},	//p4
	{2,3}	//p5
};

float oldX, oldY, oldZ, currentX, currentY, currentZ;//for bezier curve
std::vector<vector<vec3>> pointsV;//bezier curve points

std::vector<face> faceV;
std::vector<face> rotatedFaceV;
std::vector<face> transformedFaceV;

vec3 upVector = {0,1,0};

GLfloat u = 0;
GLfloat r = 4;
GLfloat alfaCam = 0;

vec3 K;	//camera location
vec3 H = {0,0,0}; //camera points towards

//uj koordrendszer
vec3 xNew, yNew, zNew;

mat4 L; //coordinate transformation matrix

GLfloat uInc = 0.3;
GLfloat rInc = 0.3;
GLfloat alfaCamInc = pi() / 45;

//light
vec3 lightSource = {1,0,0};
vec3 finalLightSource;

int factorial(int n) {
	int r = 1;
	for (int i = n; i>0; i--) {
		r *= i;
	}
	return r;
}
void drawLine(float l, float m, float j, float i)
{
	glBegin(GL_LINES);
	glVertex2f(l, m);
	glVertex2f(j, i);
	glEnd();
	glFlush();
}

float mix(float a, float b, float t)
{
	return a * (1.0f - t) + b * t;
}

float BezierQuadratic(float A, float B, float C, float t)
{
	// degree 2
	float AB = mix(A, B, t);
	float BC = mix(B, C, t);

	return mix(AB, BC, t);
}

float BezierCubic(float A, float B, float C, float D, float t)
{
	// degree 3
	float ABC = BezierQuadratic(A, B, C, t);
	float BCD = BezierQuadratic(B, C, D, t);

	return mix(ABC, BCD, t);
}

float BezierQuartic(float A, float B, float C, float D, float E, float t)
{
	// degree 4
	float ABCD = BezierCubic(A, B, C, D, t);
	float BCDE = BezierCubic(B, C, D, E, t);

	return mix(ABCD, BCDE, t);
}

float BezierFive(float A, float B, float C, float D, float E, float F, float t)
{
	// degree 5
	float ABCDE = BezierQuartic(A, B, C, D, E, t);
	float BCDEF = BezierQuartic(B, C, D, E, F, t);

	return mix(ABCDE, BCDEF, t);
}

void deCasteljau() {

	pointsV.clear();

	vec2 P0 = controlPoints[0];
	vec2 P1 = controlPoints[1];
	vec2 P2 = controlPoints[2];
	vec2 P3 = controlPoints[3];
	vec2 P4 = controlPoints[4];
	vec2 P5 = controlPoints[5];

	for (float j = 0; j <= 2 * pi(); j += vaseRotate) {
		
		vector<vec3> localVec;

		for (float t = 0; t <= 1.0; t += tInc) {

			currentX = BezierFive(P0.x, P1.x, P2.x, P3.x, P4.x, P5.x, t);
			currentY = BezierFive(P0.y, P1.y, P2.y, P3.y, P4.y, P5.y, t);

			vec3 vasePoints = { currentX, currentY, 0 };

			mat4 M = rotateY(j);

			vec4 pointH = ihToH(vasePoints);
			vec4 transformedPoint = M * pointH;

			if (transformedPoint.w != 0) {
				vec3 result = hToIh(transformedPoint);
				
				localVec.push_back(vec3(result.x, result.y, result.z));
			}
		}
		pointsV.push_back(localVec);
	}
}

void drawFace(std::vector<vec3> p, vec3 n) {

	vec3 nnS = { 0,0,perspectiveNumber };

	vec3 vecPointsToS = { nnS.x-p[0].x, nnS.y-p[0].y ,  nnS.z-p[0].z };

	GLfloat skalarProd = dot(n, nnS);

	vec3 normalizedNormalVec = normalize(n);
	
	if (skalarProd > 0) {
		normalizedNormalVec *= -1;
	}
	
	vec3 normalizedLightSource = normalize(finalLightSource);

	GLfloat greyLevel = (dot(normalizedNormalVec, normalizedLightSource) + 1) / 2;
	
	glColor3f(greyLevel, greyLevel, greyLevel);
	
	glBegin(GL_POLYGON);

	for (int i = 0; i < p.size(); i++) {
		glVertex2f(p[i].x, p[i].y);
	}
	glEnd();

	glColor3f(0, 0, 0);
	glBegin(GL_LINE_LOOP);

	for (int i = 0; i< p.size(); i++)
		glVertex2f(p[i].x, p[i].y);

	glEnd();

	//sun++

	mat4 M = w2v * projection[0];

	vec4 pointH = ihToH(finalLightSource);
	vec4 transformedPoint = M * pointH;

	vec3 result = hToIh(transformedPoint);

	glColor3f(1, 1, 0);

	GLdouble r = 4;

	glBegin(GL_POLYGON);
	for (GLdouble t = 0; t <= 2 * pi(); t += 0.01)
		glVertex2d(result.x + r * cos(t), result.y + r * sin(t));
	glEnd();
	
	//sun--
}

void calc() {

	faceV.clear();
	rotatedFaceV.clear();
	transformedFaceV.clear();

	for (int i = 0; i < pointsV.size() - 1; i++) {
		for (int j = 0; j < pointsV[i].size() - 1; j++) {

			vec3 facePoints[4];
			int countFacepoint = 0;
			for (int idx = i; idx <= i + 1; idx++)
				for (int idx2 = j; idx2 <= j + 1; idx2++) {

					vec3 vasePoint = pointsV[idx][idx2];

					facePoints[countFacepoint] = { vasePoint.x,vasePoint.y,vasePoint.z };
					countFacepoint++;

				}

			vec3 temp = facePoints[3];
			facePoints[3] = facePoints[2];
			facePoints[2] = temp;
			
			face currentFace;

			for (int i = 0; i < 4; i++) {
				currentFace.p.push_back(facePoints[i]);

			}

			currentFace.normal = cross(currentFace.p[1] - currentFace.p[0], currentFace.p[2] - currentFace.p[0]);

			faceV.push_back(currentFace);
		}
	}

	//triangles at the bottom
	for (int i = 0; i < pointsV.size(); i++) {
		vec3 currP = pointsV[i][0];
		vec3 nextP;
		if (i == pointsV.size() - 1) {
			nextP = pointsV[0][0];
		}
		else {
			nextP = pointsV[i + 1][0];
		}

		face currentFace;

		currentFace.p.push_back(currP);
		currentFace.p.push_back(nextP);
		currentFace.p.push_back(vec3{ 0,0,0 });

		currentFace.normal = cross(currentFace.p[1] - currentFace.p[0], currentFace.p[2] - currentFace.p[0]);

		faceV.push_back(currentFace);
	}

	K = { r*cos(alfaCam), u , -r * sin(alfaCam) };

	vec3 minusKHVector = { H.x - K.x,H.y - K.y ,H.z - K.z};
	
	minusKHVector *= -1;
	zNew = normalize(minusKHVector);

	vec3 skalar = cross(upVector, zNew);
	xNew = normalize(skalar);

	skalar = cross(zNew, xNew);
	yNew = normalize(skalar);

	L = coordinateTransform(K, xNew, yNew, zNew);

	//projection 
	for (int i = 0; i < faceV.size(); i++) {

		face currentFace;

		for (int j = 0; j < faceV[i].p.size(); j++) {

			vec3 vasePoint = { faceV[i].p[j].x,faceV[i].p[j].y,faceV[i].p[j].z };

			mat4 M = L;

			vec4 pointH = ihToH(vasePoint);
			vec4 transformedPoint = M * pointH;

			if (transformedPoint.w != 0) {
				vec3 result = hToIh(transformedPoint);

				currentFace.p.push_back(result);

			}
		}
		currentFace.normal = cross(currentFace.p[1] - currentFace.p[0], currentFace.p[2] - currentFace.p[0]);

		transformedFaceV.push_back(currentFace);
	}

	bool swapped = false;
	do {
		swapped = false;
		for (int i = 0; i < transformedFaceV.size() - 1; i++) {
			float avgZ1 = 0;
			float avgZ2 = 0;

			for (int j = 0; j < transformedFaceV[i].p.size(); j++) {
				avgZ1 += transformedFaceV[i].p[j].z;
			}
			for (int j = 0; j < transformedFaceV[i + 1].p.size(); j++) {
				avgZ2 += transformedFaceV[i + 1].p[j].z;
			}

			avgZ1 = perspectiveNumber - (avgZ1 / transformedFaceV[i].p.size());
			avgZ2 = perspectiveNumber - (avgZ2 / transformedFaceV[i + 1].p.size());


			if (avgZ1 > avgZ2) {

				face tempFace = transformedFaceV[i];
				transformedFaceV[i] = transformedFaceV[i + 1];
				transformedFaceV[i + 1] = tempFace;

				swapped = true;

			}
		}
	} while (swapped);

	std::vector<face> tempTransformedFaceV;

	for (int i = 0; i < transformedFaceV.size(); i++) {

		face currentFace;

		for (int j = 0; j < transformedFaceV[i].p.size(); j++) {

			vec3 vasePoint = { transformedFaceV[i].p[j].x,transformedFaceV[i].p[j].y,transformedFaceV[i].p[j].z };

			mat4 M = w2v * projection[0];

			vec4 pointH = ihToH(vasePoint);
			vec4 transformedPoint = M * pointH;

			if (transformedPoint.w != 0) {
				vec3 result = hToIh(transformedPoint);

				if (result.z == 0) {
					currentFace.p.push_back(result);
				}
			}
		}
		currentFace.normal = transformedFaceV[i].normal;
		tempTransformedFaceV.push_back(currentFace);
	}

	transformedFaceV = tempTransformedFaceV;
	//light++

	vec4 lightSourceHomogen = { lightSource.x, lightSource.y, lightSource.z, 0 };

	vec4 transformedLightSorce = transpose(inverse(L)) * lightSourceHomogen;

	finalLightSource = { transformedLightSorce.x,transformedLightSorce.y,transformedLightSorce.z};

	//light--
	

	for (int i = 0; i < transformedFaceV.size(); i++) {
		bool first = i == 1;
		drawFace(transformedFaceV[i].p, transformedFaceV[i].normal);

	}

}

void initMatrices(){
    
	scaled = scale(vec3(1.0, 1.0, 1.0));

	vec2 windowSize = { 10,10};
	vec2 windowPosition = { -5,-5 };


    vec2 viewportSize = { 450, 450 };
    vec2 viewportPosition = { winWidth / 2 - viewportSize.x / 2, winHeight / 2 - viewportSize.y / 2 };
    
	w2v = windowToViewport3(windowPosition, windowSize, viewportPosition, viewportSize);
    
	projection[0] = perspective(perspectiveNumber);
    projection[1] = ortho();
    
}

void init() {
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(0.0, winWidth, 0.0, winHeight);
    glShadeModel(GL_FLAT);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(5.0);
    glLineWidth(1.0);
    
    initMatrices();
}


void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    
    glColor3f(0, 0, 0);

	deCasteljau();
	calc();

    glutSwapBuffers();
}

void update(int value)
{
	
	if (alfaCam >= 2 * pi() || alfaCam <= -2*pi()) {
		alfaCam = 0;
	}

    glutPostRedisplay();
    
    glutTimerFunc(10, update, 0);
}


void keyPressed(unsigned char key, int x, int y)
{
	keyStates[key] = 1;


	if (keyStates['0']) {
		u -= uInc;
	}
	if (keyStates['1']) {
		u += uInc;
	}

	if (keyStates['2']) {
		r -= rInc;
	}
	if (keyStates['3']) {
		r += rInc;
	}

	if (keyStates['4']) {
		alfaCam -= alfaCamInc;
	}
	if (keyStates['5']) {
		alfaCam += alfaCamInc;
	}

	if (keyStates['6']) {
		perspectiveNumber--;
		projection[0] = perspective(perspectiveNumber);
	}
	if (keyStates['7']) {
		perspectiveNumber++;
		projection[0] = perspective(perspectiveNumber);
	}

	if (keyStates['q']) {
		lightSource.x++;
		printf("lightSource.x: %g\n", lightSource.x);
	}
	if (keyStates['w']) {
		lightSource.x--;
		printf("lightSource.x: %g\n", lightSource.x);
	}
	if (keyStates['a']) {
		lightSource.y++;
		printf("lightSource.y: %g\n", lightSource.y);
	}
	if (keyStates['s']) {
		lightSource.y--;
		printf("lightSource.y: %g\n", lightSource.y);
	}
	if (keyStates['y']) {
		lightSource.z++;
		printf("lightSource.z: %g\n", lightSource.z);
	}
	if (keyStates['x']) {
		lightSource.z--;
		printf("lightSource.z: %g\n", lightSource.z);
	}

}

void keyUp(unsigned char key, int x, int y)
{
	keyStates[key] = 0;
}

void keyOperations(int value)
{
	glutPostRedisplay();
	glutTimerFunc(10, keyOperations, 0);
}


int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(winWidth, winHeight);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("vase");
    
    init();
    glutDisplayFunc(display);

    glutTimerFunc(10, update, 0);
    
	glutTimerFunc(0, keyOperations, 0);
	glutKeyboardFunc(keyPressed);
	glutKeyboardUpFunc(keyUp);

    glutMainLoop();
    return 0;
}