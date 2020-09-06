#include <iostream>
#include <Fluid.hpp>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

#define LIMIT_VELOCITY 100.0f
#define WINDOW_NAME "DisplayFluid"

int N = 64;
int SCALE=10;
int placex = N/2;
int placy = N/2;

int dt = 1;
int T = 1000;

FluidCube* mouse_callback_cube_reference = nullptr;

void mouse_callback(int  event, int  x, int  y, int  flag, void *param)
{
	float randX = (float(rand())/float((RAND_MAX)) * LIMIT_VELOCITY) - LIMIT_VELOCITY/2;
	float randY = (float(rand())/float((RAND_MAX)) * LIMIT_VELOCITY) - LIMIT_VELOCITY/2;

	x /= SCALE;
	y /= SCALE;
	y = N - y;

	if(mouse_callback_cube_reference != nullptr)
		mouse_callback_cube_reference->addVelocity(x,y,randX,randY);
		mouse_callback_cube_reference->addDensity(x,y,100);

    if (event == EVENT_MOUSEMOVE) {
        cout << "(" << x << ", " << y << ")" << endl;
    }
}

int main() {
	// Constant results
	srand((unsigned int)time(NULL));

	namedWindow( WINDOW_NAME, WINDOW_AUTOSIZE );
	setMouseCallback(WINDOW_NAME, mouse_callback);
	Mat img1(N,N,CV_8UC3,cv::Scalar(0,0,0));
	Mat img2(N*SCALE,N*SCALE,CV_8UC3,cv::Scalar(0,0,0));


	FluidCube cube(N,0,0,dt);
	mouse_callback_cube_reference = &cube;

	for(int t=0;t<T;t++) {

		// std::cout << "Speed added" << std::endl;
		// std::cout << randX << " "<< randY << std::endl;
		// std::cout << "Velocity" << std::endl;
		// for(int i=N-1;i>=0;i--) {
		// 	for(int j=0;j<N;j++){
		// 		std::cout << cube.Vx[IX(i,j)] << " ";
		// 	}
		// 	std::cout << std::endl;
		// }
		cube.step();
		
		if(t % 10 == 0) {

			img1 = cube.getImage();
			resize(img1,img2,img2.size());
			imshow( WINDOW_NAME, img2 );
			waitKey(1);
			usleep(1000);
		}

	}

	imwrite("./abc.jpg",img2);

	return 0;
}