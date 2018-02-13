/*
ETH ZÜRICH
HPCSE II : FS2016 : VORTEX PROJECT
Authors: Franziska Krummenacher (11-917-580) & Ramona Hohl (13-816-004)
*/


#include <iostream>
#include <cstring>
#include <limits>

#include "Vortexsimulation.hpp"


int main()
{

    char inputfile[] = "./InitialConditions/vortex_C01xS5xN300xT02.bin";
    std::string f= "vortex_V01C01xS5xN300xT02.csv.";

    Vortexsimulation simulation(inputfile);
    simulation.getMesh();

    int number_of_runs = 100;
    std::cout << "starting simulation" << std::endl;
    int i=0;
    int const outputfreq = 1;

    while( i < number_of_runs && simulation.getTime() <= 0.25)
    {
    std::cout << "=======================  running run " << i << "and current voritcity time : "<<simulation.getTime()<<std::endl;
         
        if(i%outputfreq==0 || i==0)
        {
            simulation.advanceandoutput(f);
            ++i;
        }
        else
        {
          //  std::cout << "=======================  running run " << i << "and current voritcity time : "<<simulation.getTime()<<std::endl;
            simulation.advance();
            ++i;
        }

    }
    std::cout<<" run numbers : " << i <<" . for current vorticity time : "<<simulation.getTime() <<std::endl;

    std::cout<<"woohooo"<<std::endl;

    return 0;
}
