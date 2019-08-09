#ifndef TIMER
#define TIMER

#include <ctime>

/*
Simple class to calculate program timings.

Usage
-----

* call the constructor at the start of the function you want to time

* to get the current runtime between timing()-calls call timing()

* The total runtime can be accessed by calling total_timing()

*/

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
