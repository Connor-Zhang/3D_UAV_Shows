/*
Author: Yuntian Zhang
Class: ECE 6122
Date last modified: 12/2/2019

Description:

The program is a simulation of Gatech Buzzy Bowl half show simulation
using OpenGL and MPI.
Each UAV will remain on the footbal field for the first five seconds,
then flying to the (0,0,50) and maintain flying in random paths along
the surface of a virtual sphere of r=10m.

*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "iomanip"
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <GL/glut.h>
#include <chrono>
#include <thread>
#include"ECE_Bitmap.h"
using namespace std;



const int numElements = 6; // x, y, z, vx, vy, vz
double mass = 1;
const int rcvSize = 16 * 6;
double* rcvbuffer = new double[rcvSize];
double sendBuffer[numElements];

/*This is the hook Low parameter
After testing, I find when k=-1, the UAVs runs the best
*/
double k = -1 * 0.1;

/*
The following are some basic parameters of each parameter
position,velocity,force we give it and their accerlation
In order to calucalte easily, we create a variable in each direction
*/
double xPosition;
double yPosition;
double zPosition;
double xVelocity;
double yVelocity;
double zVelocity;
double xForce;
double yForce;
double zForce;
double xAccerlation;
double yAccerlation;
/*
This is the condition of each UAV during the whole process,
I will explain more about condition in the calculateuavCondition() and calculateForce() function
*/
int uavCondition = 0;
/*
Notes: in order to make it use to calculte, the zAccerlation here doesn't count the gravity.
So when calculating the force, a +10N will be counted there.
*/
double zAccerlation;

//These are the parameters related to drawing bmp file and UAV initial position
double Width = (53 + (double)1 / 3) * 0.9144;
double Length = 0.9144 * 100;
double gapX = 0.25 * Length;
double gapY = 0.5 * Width;
double LengthDraw = 0.9144 * 126;
double WidthDraw = (53 + (double)1 / 3 + 6) * 0.9144;


//Used in drawUAV function to decide the color and changing color function
double colorRed = 255;
int colorChangeCondition = -1;
int timestamp = 0; //this is for calcultating color


//Here we define the time function
void timerFunction(int id);

//Here we set for texture mapping
GLuint texture[2];
struct Image
{
	unsigned long sizeX;
	unsigned long sizeY;
	char* data;
};
typedef struct Image Image;
BMP inBitmap;
#define checkImageWidth 64
#define checkImageHeight 64
GLubyte checkImage[checkImageWidth][checkImageHeight][3];

void makeCheckImage(void) {
	int i, j, c;
	for (i = 0; i < checkImageWidth; i++) {
		for (j = 0; j < checkImageHeight; j++) {
			c = ((((i & 0x8) == 0) ^ ((j & 0x8) == 0))) * 255;
			checkImage[i][j][0] = (GLubyte)c;
			checkImage[i][j][1] = (GLubyte)c;
			checkImage[i][j][2] = (GLubyte)c;
		}
	}
}


/*
the function of drawing football field
*/
void displayFootballField()
{
	glPushMatrix();
	glBindTexture(GL_TEXTURE_2D, texture[1]);
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glBegin(GL_QUADS);
	glTexCoord2f(1, 0);
	glVertex3d(-1 * LengthDraw * 0.5, WidthDraw * 0.5, 0);
	glTexCoord2f(0, 0);
	glVertex3d(1 * LengthDraw * 0.5, WidthDraw * 0.5, 0);
	glTexCoord2f(0, 1);
	glVertex3d(1 * LengthDraw * 0.5, -1 * WidthDraw * 0.5, 0);
	glTexCoord2f(1, 1);
	glVertex3d(-1 * LengthDraw * 0.5, -1 * WidthDraw * 0.5, 0);
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
	glPopMatrix();
}


/*
this is the function of drawing Virtual Sphere, I commented in the uploaded version
But during test and debug, I uncommented it.
TAs can uncomment this to see how the UAVs flys
*/
void displayVirtualSphere()
{
	glColor3ub(255, 255, 255);
	glPushMatrix();
	glTranslatef(0, 0, 50);
	glutWireSphere(10, 50, 50);
	glPopMatrix();
}

//Here we draw the UAV based on the rank
/*
Notes: I found that the instruction ask me to change UAV color every 20 seconds, but in Piazza, the prefessor says every seconds.
So, I put two functions here are commented the first one. If TA need to see how UAV change color each 20 second, just do some change in this functions
*/
void drawUAVs()
{
	glColor3ub(colorRed, 0, 0);
	for (int i = 1; i <= 15; i++)
	{
		glPushMatrix();
		glTranslatef(rcvbuffer[i * 6], rcvbuffer[i * 6 + 1], rcvbuffer[i * 6 + 2]);
		glScalef(0.5, 0.5, 0.5);
		glTranslatef(0, 0, 0.5);
		glutSolidOctahedron();
		glPopMatrix();
	}
	
	timestamp++;
	
	/*this is the function that the color changes every time stamp*/
	if (colorRed == 255)
		colorChangeCondition = -1;
	else if (colorRed == 128)
		colorChangeCondition = 1;
	colorRed = colorRed + colorChangeCondition;

	/*this is the function that the color changes every 20 time stamps*/
	//colorChangeCondition = -1;
	//if (timestamp % 20 == 19)
	//{
	//	colorRed = colorRed + colorChangeCondition;
	//	cout<<"changed"<<endl;
	//}

}

// Reshape callback
// Window size has been set/changed to w by h pixels. Set the camera
// perspective to 45 degree vertical field of view, a window aspect
// ratio of w/h, a near clipping plane at depth 1, and a far clipping
// plane at depth 100. The viewport is the entire window.

void changeSize(int w, int h)
{
	float ratio = ((float)w) / ((float)h); // window aspect ratio
	glMatrixMode(GL_PROJECTION); // projection matrix is active
	glLoadIdentity(); // reset the projection
	gluPerspective(60.0, ratio, 0.1, 1000.0); // perspective transformation
	glMatrixMode(GL_MODELVIEW); // return to modelview mode
	glViewport(0, 0, w, h); // set viewport (drawing area) to entire window
}
/*
This is the function of displaying using openGL
*/
void renderScene()
{
	// Clear color and depth buffers
	//in order to make the scene easy to view, I set this background color
	glClearColor(0.25, 0.42, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// Reset transformations
	glLoadIdentity();
	gluLookAt(0.0, -80.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	displayFootballField();
	
	/*
	this is the function of drawing Virtual Sphere, I commented in the uploaded version
	But during test and debug, I uncommented it.
	TAs can uncomment this to see how the UAVs flys
	*/
	//displayVirtualSphere();
	
	//draw the UAVs.
	drawUAVs();
	glutSwapBuffers(); // Make it all visible
	MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
}
// mainOpenGL  - standard GLUT initializations and callbacks
void timer(int id)
{
	glutPostRedisplay();
	glutTimerFunc(100, timer, 0);
}
// mainOpenGL  - standard GLUT initializations and callbacks
void mainOpenGL(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(400, 400);

	// Send location and velocity vector in each direction
	const int numElements = 6; // x, y, z, vx, vy, vz
	const int rcvSize = 16 * 6; // (Main task + 15 UAVs) * numElements
	double* rcvbuffer = new double[rcvSize];
	double sendBuffer[numElements];

	glutCreateWindow(argv[0]);
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
	glDepthFunc(GL_LESS);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	inBitmap.read("AmFBfield.bmp");
	makeCheckImage();
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(2, texture);
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, inBitmap.bmp_info_header.width, inBitmap.bmp_info_header.height, 0,GL_BGR_EXT, GL_UNSIGNED_BYTE, &inBitmap.data[0]);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glBindTexture(GL_TEXTURE_2D, texture[1]);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, checkImageWidth, checkImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, &checkImage[0][0][0]);
	glEnable(GL_TEXTURE_2D);
	glutReshapeFunc(changeSize);
	glutDisplayFunc(renderScene);
	glutTimerFunc(100, timerFunction, 0);
	glutMainLoop();
}

// timerFunction  - called whenever the timer fires
void timerFunction(int id)
{
	glutPostRedisplay();
	glutTimerFunc(100, timerFunction, 0);
}

/*
the function of calculating normal vector
since many process like calculating force and direction require normalize vector
I make it as a function
*/
double* calculatateNormalVector(int rank)
{
	double* normalVector = new double[3];
	double xPositionTemporary, yPositionTemporary, zPositionTemporary;
	xPositionTemporary = 0 - xPosition;
	yPositionTemporary = 0 - yPosition;
	zPositionTemporary = 50 - zPosition;
	normalVector[0] = xPositionTemporary / (sqrt(xPositionTemporary * xPositionTemporary + yPositionTemporary * yPositionTemporary + zPositionTemporary * zPositionTemporary));
	normalVector[1] = yPositionTemporary / (sqrt(xPositionTemporary * xPositionTemporary + yPositionTemporary * yPositionTemporary + zPositionTemporary * zPositionTemporary));
	normalVector[2] = zPositionTemporary / (sqrt(xPositionTemporary * xPositionTemporary + yPositionTemporary * yPositionTemporary + zPositionTemporary * zPositionTemporary));
	return normalVector;
}

/*
this is the function of calcultating the UAV's current position to the center
of the sphere.
I put it as a function since lots of place reqirues this kind of calculation
*/
double calculateDistance(int rank)
{
	return sqrt((xPosition - 0) * (xPosition - 0) + (yPosition - 0) * (yPosition - 0) + (zPosition - 50) * (zPosition - 50));
}

/*
This function is calucating the Uav's position and velocity function
After assigning force (in each direction) new position and velocity is calculated via function
the new information are sent to sendBuffer for next round's gathering inforation to the main thread
*/
void CalcualteUAVsInformation(int rank)
{
	xPosition = xPosition + xVelocity * 0.1 + 0.5 * xAccerlation * 0.01;
	yPosition = yPosition + yVelocity * 0.1 + 0.5 * yAccerlation * 0.01;
	zPosition = zPosition + zVelocity * 0.1 + 0.5 * zAccerlation * 0.01;
	xVelocity = xVelocity + xAccerlation * 0.1;
	yVelocity = yVelocity + yAccerlation * 0.1;
	zVelocity = zVelocity + zAccerlation * 0.1;
	sendBuffer[0] = xPosition;
	sendBuffer[1] = yPosition;
	sendBuffer[2] = zPosition;
	sendBuffer[3] = xVelocity;
	sendBuffer[4] = yVelocity;
	sendBuffer[5] = zVelocity;
}

/*
in the program, I design seven different conditions for the UAVs.
this function is used to each condition's function
I will explain more about each condition and their change conditions in the calculate Force() function
*/
void calculateUavCondition(int rank)
{
	double tempDistance = calculateDistance(rank);
	double* normalVector = calculatateNormalVector(rank);
	double tempCurrentlyVelocity = sqrt(xVelocity*xVelocity + yVelocity*yVelocity + zVelocity*zVelocity);
	if (uavCondition == 0)
	{
		if (tempCurrentlyVelocity > 1.9)
			uavCondition = 1;
	}
	else if (uavCondition == 1)
	{
		if ((tempDistance - 2) < 10)
			uavCondition = 2;
	}
	else if (uavCondition == 2)
	{
		if (tempCurrentlyVelocity < 0.1)
			uavCondition = 3;
	}
	else if (uavCondition == 4)
	{
		if (tempCurrentlyVelocity > 2.4)
		{
			uavCondition = 5;
		}
	}
	else if (uavCondition == 5)
	{
		if (tempCurrentlyVelocity > 2.6)
		{
			uavCondition = 6;
		}
		if (tempCurrentlyVelocity < 2.4)
		{
			uavCondition = 6;
		}
	}
	else if (uavCondition == 6)
	{
		if (tempCurrentlyVelocity < 2.6)
		{
			uavCondition = 5;
		}
		if (tempCurrentlyVelocity > 2.4)
		{
			uavCondition = 5;
		}
	}
}
/*
After deciding the Uav condition, this funiction is used to assign force to each
Uav based on their current condition.
There are six conditions in my program, I will explain it in the function
of how each conditions mean and there how to assign force to them
*/
void calculateForce(int rank)
{
	double tempDistance = calculateDistance(rank);
	double* normalVector = calculatateNormalVector(rank);
	double tempCurrentlyVelocity = sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity);
	double tempHookForce = k * (10 - tempDistance) + 0.625;
	/*
	Condition 0 is the first condition of each UAV each are assgined this at first
	In this condition, the UAV are trying to go to the sphere in a accelerating speed
	Therefore, I give each direction's force twice as the normalVector to give them force to speed up
	Notes:To make sure the maximum force is less then 20N, I assign the multiplier as 2 here.
	*/
	if (uavCondition == 0)
	{
		xForce = 2*normalVector[0];
		yForce = 2*normalVector[1];
		zForce = 2*normalVector[2] + 10;
		xAccerlation = xForce / mass;
		yAccerlation = yForce / mass;
		zAccerlation = (zForce - 10) / mass;
	}
	/*
	Condition 1 is when we detecting the velocity of Uav has reached to some point during its flying to
	the sphere, we don't put any other force (except the force to eliminiate gravity).
	In this condition, the UAV are moving to the sphere in a constant spped.
	*/
	else if (uavCondition == 1)
	{
		xForce = 0;
		yForce = 0;
		zForce = 0;
		xAccerlation = 0;
		yAccerlation = 0;
		zAccerlation = 0;
	}
	/*
	when we detect the UAV is about to land the surface of the sphere (in a current distance),
	we move uavCondition into condition 2.
	This is a condition that the Uav is appoaching sphere in a decelerating speed.
	Opposite as Condition 0 ,we give a negative force along the uav's moving direction to make it slow thow.
	This step is designed for steps.
	*/
	else if (uavCondition == 2)
	{
		xForce = -2 * normalVector[0];
		yForce = -2 * normalVector[1];
		zForce = -2 * normalVector[2] +10;
		xAccerlation = xForce / mass;
		yAccerlation = yForce / mass;
		zAccerlation = (zForce - 10) / mass;
	}
	/*
	This is a condition at the moment that the Uav reaches the surface of the sphere
	we make a random seed to assgin a tangential force to the Uav
	A normal vector is counted first to do this step
	Since this condition is temporary and immediate, we assign the condition to condition 4 at last
	and don't do any condition logic in the CalculateUavCondition() function
	*/
	else if (uavCondition == 3)
	{
		srand((unsigned int)time(NULL));
		double xTemporyTangentialVector = 2*((double)rand() / RAND_MAX) -1;
		double yTemporyTangentialVector = 2*((double)rand() / RAND_MAX)-1;
		double zTemporyTangentialVector = -1 * (xTemporyTangentialVector * normalVector[0] +yTemporyTangentialVector * normalVector[1]) / normalVector[2];
		double xTangentialVector = xTemporyTangentialVector / (sqrt(xTemporyTangentialVector *xTemporyTangentialVector +yTemporyTangentialVector * yTemporyTangentialVector + zTemporyTangentialVector * zTemporyTangentialVector));
		double yTangentialVector = yTemporyTangentialVector / (sqrt(xTemporyTangentialVector *xTemporyTangentialVector +yTemporyTangentialVector * yTemporyTangentialVector + zTemporyTangentialVector * zTemporyTangentialVector));
		double zTangentialVector = zTemporyTangentialVector / (sqrt(xTemporyTangentialVector *xTemporyTangentialVector +yTemporyTangentialVector * yTemporyTangentialVector + zTemporyTangentialVector * zTemporyTangentialVector));
		double normalDeviationForce = (tempCurrentlyVelocity * tempCurrentlyVelocity) / tempDistance;
		xForce = xTangentialVector + (normalDeviationForce + k * (10 - tempDistance)) * normalVector[0];
		yForce = yTangentialVector + (normalDeviationForce + k * (10 - tempDistance)) * normalVector[1];
		zForce = zTangentialVector + (normalDeviationForce + k * (10 - tempDistance)) * normalVector[2] + 10;
		xAccerlation = xForce / mass;
		yAccerlation = yForce / mass;
		zAccerlation = (zForce - 10) / mass;
		/*
		Since condition 3 is a immediate process, we put the condition as uavCondition4 immediately after the force is assgined
		*/
		uavCondition = 4;
	}
	/*
	the Uav is immeaditely moved into condition 4 when they are assigned a tangential force at condition 3
    The ideal speed of UAV is between 2.4-2.6. So in Condition4, our force given to each Uav is a little bit higher than
	it should be when moving in a constant speed: to make it speed up.
	When the speed is higher than 2.4, uav will come to condition 5
	*/
	else if (uavCondition == 4)
	{
		double normalDeviationForce = (tempCurrentlyVelocity * tempCurrentlyVelocity) / tempDistance;
		xForce = xVelocity / (sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity)) + (normalDeviationForce + k * (10 - tempDistance)) * normalVector[0];
		yForce = yVelocity / (sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity)) + (normalDeviationForce + k * (10 - tempDistance)) * normalVector[1];
		zForce = zVelocity / (sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity)) + (normalDeviationForce + k * (10 - tempDistance)) * normalVector[2] + 10;
		xAccerlation = xForce / mass;
		yAccerlation = yForce / mass;
		zAccerlation = (zForce - 10) / mass;
	}
	/*
	Condition 5 is the most ideal condition: that the uav is moving in a nearly constant speed that is approximate
	to what I want the uav is moving. 
	So in this process, the force given is just the hook distance and the gravity
	*/
	else if (uavCondition == 5)
	{
		double dir = -1;
		if (abs(tempHookForce) < 9)
		{
			xForce = normalVector[0] * (k * (10 - tempDistance) + 0.625);
			yForce = normalVector[1] * (k * (10 - tempDistance) + 0.625);
			zForce = normalVector[2] * (k * (10 - tempDistance) + 0.625) + 10;
		}
		else
		{
			if (tempHookForce > 0)
			{
				dir = 1;
			}
			xForce = normalVector[0] * (9 * dir + 0.625);
			yForce = normalVector[1] * (9 * dir + 0.625);
			zForce = normalVector[2] * (9 * dir + 0.625) + 10;
		}
		xAccerlation = xForce / mass;
		yAccerlation = yForce / mass;
		zAccerlation = (zForce - 10) / mass;
	}
	/*
	Condition 6 is a condition that the uav is a little bit out of the reasonable speed to let the 
	uav move constantly(2.4-2.6) 
	When we detect the UAV is lower or higher than the range of reasonable speed, I move them to condition 6
	In condtion 6, more force are given(compare to condition 5) to make up/offset the speed more than/less than
	the reasonable speed either positively or negatively.
	*/
	else if (uavCondition == 6)
	{
		double dir;
		if (tempCurrentlyVelocity > 2.6)
			dir = -1;
		else if (tempCurrentlyVelocity < 2.4)
			dir = 1;
		else
			dir = 0;
		xForce = xVelocity / (sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity)) *0.4 * dir + normalVector[0] * 0.625;
		yForce = yVelocity / (sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity)) *0.4 * dir + normalVector[1] * 0.625;
		zForce = zVelocity / (sqrt(xVelocity * xVelocity + yVelocity * yVelocity + zVelocity * zVelocity)) * 0.4 * dir + normalVector[2] * 0.625 + 10;
		xAccerlation = xForce / mass;
		yAccerlation = yForce / mass;
		zAccerlation = (zForce - 10) / mass;
	}
}




// Main entry point determines rank of the process and follows the correct program path
int main(int argc, char** argv)
{
	int numTasks, rank;
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS)
	{
		printf("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int gsize = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &gsize);

	
	/*
	this is a 2D vector storing the initial information (position & velocity of each UAVs.
	For convience, I assigned 16 vectors here to match the number of vector to the rank number
	so the first vector is useless is just for test use.
	*/
	double uavInitialInformation[16][6] = { { 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 }, //test no use
	{ -2 * gapX, gapY, 0.00, 0.0, 0.0, 0.0 }, { -1 * gapX, gapY, 0.0, 0.0, 0.0, 0.0 }, { 0.0, gapY, 0.0, 0.0, 0.0, 0.0 }, { 1 * gapX, gapY, 0.0, 0.0, 0.0, 0.0 }, { 2 * gapX, gapY, 0.0, 0.0, 0.0, 0.0 },
	{ -2 * gapX, 0, 0.00, 0.0, 0.0, 0.0 }, { -1 * gapX, 0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0, 0.0, 0.0, 0.0, 0.0 }, { 1 * gapX, 0, 0.0, 0.0, 0.0, 0.0 }, { 2 * gapX, 0, 0.0, 0.0, 0.0, 0.0 },
	{ -2 * gapX, -gapY, 0.00, 0.0, 0.0, 0.0 }, { -1 * gapX, -gapY, 0.0, 0.0, 0.0, 0.0 }, { 0.0, -gapY, 0.0, 0.0, 0.0, 0.0 }, { 1 * gapX, -gapY, 0.0, 0.0, 0.0, 0.0 }, { 2 * gapX, -gapY, 0.0, 0.0, 0.0, 0.0 },
	};
	
	if (rank == 0)
	{
		for (int i = 0; i < 6; i++)
		{
			sendBuffer[i] = uavInitialInformation[rank][i];
		}
		xAccerlation = 0.0;
		yAccerlation = 0.0;
		zAccerlation = 0.0;
		xForce = 0.0;
		yForce = 0.0;
		zForce = 0.0;
		xPosition = sendBuffer[0];
		yPosition = sendBuffer[1];
		zPosition = sendBuffer[2];
		xVelocity = sendBuffer[3];
		yVelocity = sendBuffer[4];
		zVelocity = sendBuffer[5];
	}
	else
	{
		for (int i = 0; i < 6; i++)
		{
			sendBuffer[i] = uavInitialInformation[rank][i];
		}
		for (int i = 0; i < rcvSize; i++)
		{
			rcvbuffer[i] = 0;
		}
		xPosition = sendBuffer[0];
		yPosition = sendBuffer[1];
		zPosition = sendBuffer[2];
		xVelocity = sendBuffer[3];
		yVelocity = sendBuffer[4];
		zVelocity = sendBuffer[5];
		xAccerlation = 0.0;
		yAccerlation = 0.0;
		zAccerlation = 0.0;
		xForce = 0.0;
		yForce = 0.0;
		zForce = 0.0;
	}

	/*
	After assigning each rank its Uav's initial information
	The ALlgather is used to gather and broadcase all Uav information
	*/
	MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);

	//Here the simulation begins
	if (rank == 0)
	{
		mainOpenGL(argc, argv);
	}
	else
	{
		// Sleep for 5 seconds
		std::this_thread::sleep_for(std::chrono::seconds(5));
		/*
		For each timestamp, I first calculate the Uav condition, then assigning them force in 
		each direction based on their conditions. Finally update their position and velocity information
		and gather & broadcase to all others.
		*/
		/*
		In my test, the first sphere will get to the surfect of sphere at about 232 timestamp, so 
		I set the whole timestamp as 832 timestamp
		*/
		for (int ii = 0; ii < 832; ii++)
		{

			calculateUavCondition(rank);
			calculateForce(rank);
			CalcualteUAVsInformation(rank);

			MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE,rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
			
			/*
			this is the function of testing if the two UAVs will bump to each other and collise
			if so, we sawp their velocity
			*/
			for (int i = 1; i < 16; i++)
			{
				if (i == rank)
					continue;
				if (sqrt((xPosition - rcvbuffer[i * 6]) * (xPosition - rcvbuffer[i * 6])
					+ (yPosition - rcvbuffer[i * 6 + 1])*(yPosition - rcvbuffer[i * 6 + 1])
					+ (zPosition - rcvbuffer[i * 6 + 2])*(zPosition - rcvbuffer[i * 6 + 2]))<1.01
					&& 
					uavCondition > 2)
				{
					xVelocity = rcvbuffer[i * 6 + 3];
					yVelocity = rcvbuffer[i * 6 + 4];
					zVelocity = rcvbuffer[i * 6 + 5];
				}
			}
		}
	}
	return 0;
}
