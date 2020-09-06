#include <iostream>
#include <Fluid.hpp>
#include <cstdlib>
#include <ctime>

#define LIMIT_VELOCITY 100.0f

int main() {
	// Constant results
	srand((unsigned int)time(NULL));

	int N = 64;
	int dt = 0.1;

	FluidCube cube(N,0,0,dt);
	cube.addDensity(2,2,100);

	for(int t=0;t<1000;t++) {

		float randX = (float(rand())/float((RAND_MAX)) * LIMIT_VELOCITY);
		float randY = (float(rand())/float((RAND_MAX)) * LIMIT_VELOCITY);

		cube.addVelocity(20,20,randX,randY);

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
	}

	cube.saveImage("./abc.jpg");

	return 0;
}