// -*- C++ -*-
//	file:	wire2d.cpp
//	brief:  Simple realtime 2D wire solver
//	usage:  wire2d [-l <float>] [-n <int>] [-s <float>] [-g <float>]
//		-l <float>(default 200)		wire length in pixels.
//	 	-n <int>(default 10)		number of wire segments
//		-s <float>(default 1)		wire stiffness. Possible values are in range(0,1]
//		-g <float>(default 20)		gravity force
//
// 
// Copyright (c) 2016 Roman Saldygashev <sldg.roman@gmail.com>
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#include <SDL2/SDL.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cmath>


// helper function for logging
void
log(const std::string && message)
{
    std::cout << message << std::endl;
}

// helper function for error logging
void
error_log(const std::string && message)
{
    std::cerr << message << std::endl
	      << "SDL Error:\n"
	      << SDL_GetError() << std::endl;
}


// generic two dimensional vector
template <typename T>
struct Vector2
{
    T x;
    T y;

    Vector2
    operator+(const Vector2 & rhs) const noexcept
    { return Vector2{x+rhs.x, y+rhs.y}; }
    Vector2
    operator+(const T & rhs) const noexcept
    { return Vector2{x+rhs, y+rhs}; }
    Vector2&
    operator+=(const Vector2 & rhs)
    { *this = *this + rhs; return *this; }
    Vector2&
    operator+=(const T & rhs)
    { *this = *this + rhs; return *this; }

    Vector2
    operator-(const Vector2 & rhs) const noexcept
    { return Vector2{x-rhs.x, y-rhs.y}; }
    Vector2
    operator-(const T & rhs) const noexcept
    { return Vector2{x-rhs, y-rhs}; }
    Vector2
    operator-() const noexcept
    { return Vector2{-x, -y}; }
    Vector2&
    operator-=(const Vector2 & rhs)
    { *this = *this - rhs; return *this; }
    Vector2&
    operator-=(const T & rhs)
    { *this = *this - rhs; return *this; }

    Vector2
    operator*(const Vector2 & rhs) const noexcept
    { return Vector2{x*rhs.x, y*rhs.y}; }
    Vector2
    operator*(const T & rhs) const noexcept
    { return Vector2{x*rhs, y*rhs}; }
    Vector2&
    operator*=(const Vector2 & rhs)
    { *this = *this * rhs; return *this; }
    Vector2&
    operator*=(const T & rhs)
    { *this = *this * rhs; return *this; }

    Vector2
    operator/(const Vector2 & rhs) const noexcept
    { return Vector2{x/rhs.x, y/rhs.y}; }
    Vector2
    operator/(const T & rhs) const noexcept
    { return Vector2{x/rhs, y/rhs}; }

    template <typename T2>
    operator Vector2<T2>() const noexcept
    { return Vector2<T2>{T2(x), T2(y)}; }
};


// geometry functions for Vector2
template <typename T>
T
length(const Vector2<T> & a)
{ return std::sqrt(a.x*a.x+a.y*a.y); }

template <typename T>
Vector2<T>
normalize(const Vector2<T> & a)
{ return a / length(a); }

template <typename T>
T
dot(const Vector2<T> & lhs, const Vector2<T> & rhs)
{ return lhs.x*rhs.x + lhs.y*rhs.x; }

template <typename T>
T
clamp(T val, const T & l, const T & h)
{
    val = val < l ? l : val;
    return val > h ? h : val;
}


using Vec2f = Vector2<float>;
using Vec2i = Vector2<int>;


// class for the time increment calculation
class Timer
{
    int prev_time = SDL_GetTicks();

public:
    Timer() {}

    float
    get_time_increment()
    {
	int ret = SDL_GetTicks() - prev_time;
	prev_time = ret;
	return ret;	// convert milliseconds to seconds
    }

    int leftOver = 0;
};


const float minf = 1.5f;
const float maxf = 10.f;


// 2-dimensional point
// assume the mass always is 1
struct Point
{
    Vec2f p; 			// current position
    Vec2f pp;			// previous position

    void
    move(const Vec2f & pos) noexcept
    {
    	pp = p;
    	p = pos;
    }
};


// simple link constraint
struct Link
{
    float rest; 		// rest distance
    float stiffness;		// spring stiffness

    Point* p1;
    Point* p2;

    // solve only point 2
    void
    solve_p2() noexcept
    {
	Vec2f dir = p1->p - p2->p;
	float dist = length(dir);
	float f = std::abs(rest - dist) / dist;
	if (rest >= dist)
	    dir = -dir;

	p2->pp = p2->p;
	p2->p += normalize(dir) * f * stiffness;
    }

    // solve both points
    void
    solve() noexcept
    {
	Vec2f dir = p1->p - p2->p;
	float dist = length(dir);
	float f = std::abs(rest - dist) / dist;
	if (rest >= dist)
	    dir = -dir;

	p1->pp = p1->p;
	p1->p -= normalize(dir) * f * stiffness * 0.5;
	p2->pp = p2->p;
	p2->p += normalize(dir) * f * stiffness * 0.5;
    }
};


// simple position-based solver with Verlet integration
Point&
solve_point(Point & point, const float dt, const Vec2f F, const Vec2f & ss)
{
    Vec2f v = point.p - point.pp;

    if (point.p.x < 1) {
	point.p.x = 1;
	if (v.x < 0)
	    v.x *= -1.f;
    }
    if (point.p.y < 1) {
	point.p.y = 1;
	if (v.y < 0)
	    v.y *= -1.f;
    }
    if (point.p.x > ss.x-1) {
	point.p.x = ss.x-1;
	if (v.x > ss.x-1)
	    v.x *= -1.f;
    }
    if (point.p.y > ss.y-1) {
	point.p.y = ss.y-1;
	if (v.y > ss.y)
	    v.y *= -1.f;
    }

    Vec2f newp = point.p * 2.f - point.pp + v + F * dt*dt;
    point.pp = point.p;
    point.p = newp;

    return point;
}

// default resolution
const int xres = 640;
const int yres = 480;

// wrapper for SDL
class SDLWrapper
{
    Vec2i screen_size {xres, yres};

    SDL_Window* gWindow = nullptr;
    SDL_Renderer* gRenderer = nullptr;
    bool _status = true;

public:
    SDLWrapper()
    {
	if( SDL_Init( SDL_INIT_VIDEO ) < 0 )
	    {
		error_log("SDL could not initialize!");
		_status = false;
	    }

	gWindow = SDL_CreateWindow("Wire",
				   SDL_WINDOWPOS_UNDEFINED,
				   SDL_WINDOWPOS_UNDEFINED,
				   screen_size.x,
				   screen_size.y,
				   SDL_WINDOW_SHOWN);
	if(gWindow == nullptr)
	    {
		error_log("Window could not be created!");
		_status = false;
	    }

	gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
	if (gRenderer == nullptr)
	    {
		error_log("Renderer could not be created!");
		_status = false;
	    }

	SDL_SetRenderDrawColor(gRenderer, 0x00, 0x00, 0x00, 0x00);

	log("Renderer has initialized.");
    }

    // dynamic window resize support
    void
    update_window() noexcept
    { SDL_GetWindowSize(gWindow, &(screen_size.x), &(screen_size.y)); }

    bool
    status() noexcept
    { return _status; }

    void
    set_black() noexcept
    { SDL_SetRenderDrawColor(gRenderer, 0x00, 0x00, 0x00, 0xFF); }

    void
    clear() noexcept
    { SDL_RenderClear(gRenderer); }


    void
    render() noexcept
    { SDL_RenderPresent(gRenderer); }

    void
    render_clear() noexcept
    { SDL_RenderClear(gRenderer); }

    void
    draw_line(const int x1, const int y1, const int x2, const int y2) noexcept
    {
	SDL_SetRenderDrawColor(gRenderer, 0xFF, 0xFF, 0xFF, 0xFF);
	SDL_RenderDrawLine(gRenderer, x1, y1, x2, y2);
    }

    void
    draw_point(const int x, const int y)
    {
	SDL_SetRenderDrawColor(gRenderer, 0xFF, 0x00, 0x00, 0xFF);
	SDL_RenderDrawPoint(gRenderer, x, y);
    }

    void
    draw_link(const Link & link)
    {
	draw_line(link.p1->p.x, link.p1->p.y,
		  link.p2->p.x, link.p2->p.y);
    }


    Vec2i
    get_screen_size() const noexcept
    { return screen_size; }

    virtual
    ~SDLWrapper()
    {
	SDL_DestroyRenderer(gRenderer);
	SDL_DestroyWindow(gWindow);
	gWindow = nullptr;
	gRenderer = nullptr;
	SDL_Quit();

	log("Renderer has destroyed.");
    }
};



int main(int argc, char *argv[])
{
    // default values
    float wirelength = 200.f;	// wire length in pixels
    size_t nlinks = 20; 	// number of segments
    float stiffness = 1.f; 	// wire stiffness
    float gravity = 20.f;	// gravity force

    // argument handling
    if (argc < 2)
	std::cout << "Usage: wire2d "
		  << "[-l <float>] [-n <int>] [-s <float>] [-g <float>]\n"
		  << "\t-l <float>(default " << wirelength << ")\t\twire length in pixels.\n"
		  << "\t-n <int>(default " << nlinks << ")\t\tnumber of wire segments\n"
		  << "\t-s <float>(default " << stiffness << ")\t\twire stiffness. Possible values are in range(0.1,1]\n"
		  << "\t-g <float>(default " << gravity << ")\t\tgravity force\n";
    else {
	if (argc % 2 != 1)
	    log("Wrong number of arguments!");
	else
	    for (int i = 1; i < argc-1; i+=2) {
		std::string flag(argv[i]);
		std::string val(argv[i+1]);
		if (flag == "-l") wirelength = std::abs(std::stof(val));
		else if (flag == "-n") nlinks = std::stoi(val);
		else if (flag == "-s") stiffness = clamp(std::stof(val), 0.1f, 1.f);
		else if (flag == "-g") gravity = std::stof(val);
		else std::cout << "Unknown argument: " << flag << std::endl;
	    }
    }
    // std::cout << "Length = " << wirelength << std::endl;

    // points
    size_t npoints = nlinks + 1;
    std::vector<Point> points;
    points.reserve(npoints);
    // // randomize initial points positions
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> xdstrb(1, xres);
    std::uniform_real_distribution<float> ydstrb(1, yres);
    std::generate_n(std::back_inserter(points), npoints,
		    [&](){ return Point{Vec2f{xdstrb(gen), ydstrb(gen)},
				Vec2f{xdstrb(gen), ydstrb(gen)}}; });
    points.shrink_to_fit();

    // links
    float segment_length = wirelength / nlinks;
    std::vector<Link> links;
    links.reserve(nlinks);
    for (size_t i = 0; i < nlinks; i++) {
	links.push_back(Link{segment_length, stiffness,
		    &(points[i]), &(points[i+1])});
    }
    links.shrink_to_fit();

    Point masterPoint {Vec2f{xdstrb(gen), ydstrb(gen)},
	    Vec2f{xdstrb(gen), ydstrb(gen)}};
    Link masterLink {segment_length, stiffness,
	    &masterPoint, &(points[0])};

    // gravity force
    Vec2f g {0.f, gravity};

    // system initialization
    SDLWrapper renderer;

    if (!renderer.status())
	return 1;

    bool quit = false;
    SDL_Event event;
    Timer timer;

    // game loop
    while (!quit)
	{
	    while(SDL_PollEvent(&event) != 0)
		if (event.type == SDL_QUIT)
		    quit = true;

	    Vec2i mouse_coord;

	    if (event.type == SDL_MOUSEMOTION)
	        SDL_GetMouseState(&(mouse_coord.x),
				  &(mouse_coord.y));

	    // mouse interaction
	    masterPoint.move(Vec2f(mouse_coord));

	    Vec2f ss = Vec2f(renderer.get_screen_size());

	    int elapsed = timer.get_time_increment();
	    elapsed += timer.leftOver;
	    int timesteps = std::floor(elapsed / 16);
	    timer.leftOver = elapsed - timesteps * 16;
	    float dt = float(elapsed) / float(timesteps) * 0.001;

	    // solve dynamics
	    for (int i = 0; i < timesteps; i++) {
		masterLink.solve_p2();

		for (auto & link : links) {
		    // solve points first
		    *link.p1 = solve_point(*link.p1, dt, g, ss);
		    *link.p2 = solve_point(*link.p2, dt, g, ss);
		    // solve link next
		    link.solve();
		}
	    }

	    // render result
	    renderer.set_black();
	    renderer.clear();

	    renderer.draw_link(masterLink);
	    for (auto link : links){
		renderer.draw_link(link);
	    }

	    renderer.render();
	    renderer.update_window();
	}

    return 0;
}
