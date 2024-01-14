#include <iostream>
#include "Scheduler.h"

void Scheduler::insertEvent(double lambda, const scheduleEvent &e)
{
    schedule.insert(pair<double, scheduleEvent>(lambda, e));
}

pair<double, scheduleEvent> Scheduler::getNextEvent(bool forceOrder)
{
    pair<double, scheduleEvent> foo;
    multimap<double, scheduleEvent>::iterator iter, iter2;
    
    // get an iterator for the first element
    if(forceOrder == false) {
        iter = schedule.begin();
        foo = *iter;
        schedule.erase(iter);
        return(foo);
    }
    else { // use priority ordering; merging > {tension, activate, deactivate}
        iter = schedule.begin();
        if(iter->second.type !='M') { // not a merge event; search if there is a merge event within the tolerance level
            double lambda = iter->first;
            for(iter2=iter; iter2!=schedule.end() && iter2->first < lambda+tolerance; ++iter2) {
                if(iter2->second.type == 'M') { iter=iter2;}
            }
        }
        foo = *iter;
        schedule.erase(iter);
        return(foo);
    }
}

// clears the whole scheduler of all events
void Scheduler::clearSchedule()
{
    schedule.clear();
}




void Scheduler::printSchedule(ostream& outStream)
{
    multimap<double, scheduleEvent>::iterator iter;
    // go through all scheduled events
    for(iter = schedule.begin(); iter!=schedule.end(); ++iter)
    {
        outStream << "Lambda: " << iter->first << endl;
        outStream << "Type: " << iter->second.type << " Group 1: " << iter->second.grp1 << " Group2: " << iter->second.grp2 << endl;
    }
    outStream << endl;
}
