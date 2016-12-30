#ifndef TIMER
#define TIMER

#include <ctime>


//creates object that can keep track of time in program
class timer 
{

	private:
		clock_t start_stamp, last_stamp, current_stamp;	

	public:
		timer() 
		{ 
			start_stamp = clock();
			last_stamp = start_stamp;
		}
		float timing() 
		{
			current_stamp = clock();
			float t_diff ((float)current_stamp - (float)last_stamp);
			last_stamp = current_stamp;
			return t_diff/CLOCKS_PER_SEC;
		}
		float total_timing() 
		{
			current_stamp = clock();
			float t_diff ((float)current_stamp - (float)start_stamp);
			return t_diff/CLOCKS_PER_SEC;
		}
};


#endif
