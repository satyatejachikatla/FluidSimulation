#include <iostream>
#include <vector>
#include <math.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <Fluid.hpp>

void foo() {
	std::cout << "Foo" << std::endl;
}

FluidCube::FluidCube(int size, int diffusion, int viscosity, float dt) {

	int N = size;

	this->size = size;
	this->dt = dt;
	this->diff = diffusion;
	this->visc = viscosity;

	this->s.resize(N * N * N,0);
	this->density.resize(N * N * N,0);

	this->Vx.resize(N * N * N,0);
	this->Vy.resize(N * N * N,0);

	this->Vx0.resize(N * N * N,0);
	this->Vy0.resize(N * N * N,0);

}

FluidCube::~FluidCube() {

}

void FluidCube::addDensity(int x, int y, float amount) {

    int N = this->size;
    int index = IX(x, y);

	this->density[index] += amount;
}

void FluidCube::addVelocity(int x, int y, float amountX, float amountY) {

    int N = this->size;
	int index = IX(x, y);

	this->Vx[index] += amountX;
	this->Vy[index] += amountY;

}

static void set_bnd(int b,std::vector<float>& x, int N) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1  )];
            x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
        }
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
            x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
        }

    x[IX(0, 0)]         = 0.5 * (x[IX(1, 0)]         + x[IX(0, 1)]);
    x[IX(0, N - 1)]     = 0.5 * (x[IX(1, N - 1)]     + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)]     = 0.5 * (x[IX(N - 2, 0)]     + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);

}

static void lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c, int iter, int N) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a*(    x[IX(i+1, j)]
                                    +x[IX(i-1, j)]
                                    +x[IX(i  , j+1)]
                                    +x[IX(i  , j-1)]
                           )) * cRecip;
                }
            }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, std::vector<float>& x, std::vector<float>& x0, float diff, float dt, int iter, int N) {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div, int iter, int N) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                         velocX[IX(i+1, j)]
                        -velocX[IX(i-1, j)]
                        +velocY[IX(i  , j+1)]
                        -velocY[IX(i  , j-1)]
                    )/N;
                p[IX(i, j)] = 0;
            }
        }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);
    
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N;
            }
        }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

static void advect(int b, std::vector<float>& d,std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY, float dt, int N) {
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floorf(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floorf(y);
                j1 = j0 + 1.0f;
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
                
                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;
                
                d[IX(i, j)] = 
                
                    s0 * ( t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                    s1 * ( t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
            }
        }
    set_bnd(b, d, N);
}

void FluidCube::step() {
	int N = this->size;
    int iter = 4;

	diffuse(1, Vx0, Vx, visc, dt, iter, N);
	diffuse(2, Vy0, Vy, visc, dt, iter, N);

	project(Vx0, Vy0, Vx, Vy, iter, N);

	advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
	advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

	project(Vx, Vy, Vx0, Vy0, iter, N);

	diffuse(0, s, density, diff, dt, iter, N);
	advect(0, density, s, Vx, Vy, dt, N);
}


float clamp(float x0,float x1,float n) {
    n = abs(n);
    if ( x0 <= n && n <= x1 )
        return n;
    if ( n <= x0 )
        return x0;
    return x1;
}

void FluidCube::saveImage(std::string filePath){

    int N = this->size;
    int nx,ny;
    nx=ny=N;

    cv::Mat img(ny,nx,CV_8UC3,cv::Scalar(0,0,0));

    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {

            int ir,ig,ib;
            ir  = clamp(0,255,Vx[IX(i,j)]);
            ig  = clamp(0,255,Vy[IX(i,j)]);
            ib  = clamp(0,255,density[IX(i,j)]);

            //ir = ig = ib;

            img.at<cv::Vec3b>(ny-j-1,i)[0] = ib;/*B*/
            img.at<cv::Vec3b>(ny-j-1,i)[1] = ig;/*G*/
            img.at<cv::Vec3b>(ny-j-1,i)[2] = ir;/*R*/
        }
    }

    cv::imwrite(filePath,img);
}

