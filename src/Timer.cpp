#include "../hdr/Timer.h"
#include <iostream>
#include <cmath>
#include <ctime>
using namespace std;

int main(int argc, char *argv[])
{
    double start, stop, tick1, tick2;

    // second, use Timer::getElapsedTime() ////////////////////////////////////
    Timer t;

    // start timer
    t.start();
    tick1 = tick2 = t.getElapsedTimeInMilliSec();

    while(t.getElapsedTime() < 1)       // loop for 1 sec
    {
        cout << (tick2) << " ms." << endl;

        tick1 = tick2;
        tick2 = t.getElapsedTimeInMilliSec();
    }

    cout << t.getCurrentDatetime() << endl;

    return 0;
}
