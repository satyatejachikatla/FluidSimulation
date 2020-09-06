#pragma once

#include <vector>
#include <string>

#define IX(x, y) ((x) + (y) * N)

class FluidCube {
	private:
	public:
		int size;
		float dt;
		float diff;
		float visc;

		std::vector<float> s;

		std::vector<float> Vx;
		std::vector<float> Vy;

		std::vector<float> Vx0;
		std::vector<float> Vy0;

	
		std::vector<float> density;

		FluidCube(int size, int diffusion, int viscosity, float dt);
		~FluidCube();

		void addDensity(int x, int y, float amount);
		void addVelocity(int x, int y, float amountX, float amountY);
		void step();

		void saveImage(std::string filePath);
};