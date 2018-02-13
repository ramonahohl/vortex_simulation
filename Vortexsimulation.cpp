#include "Vortexsimulation.hpp"


Vortexsimulation::Vortexsimulation(const char filename[])
{
    initialize(filename);
    dx_ = (xmax_ - xmin_)/(M_-1);
    cutoff_ = 10.*dx_;
    T_ = 0.;
    outnumber_=0;
}

Vortexsimulation::~Vortexsimulation()
{
    free(xsources_);
    free(ysources_);
    free(wsources_);
    free(xmesh_);
    free(ymesh_);
    free(phi_);
    free(ux_);
    free(uy_);
    free(vorticity_);
    free(vorticity_temp_);
}

void Vortexsimulation::initialize(const char filename[])
{

    char file[256];
    strcpy(file, filename);

    if (access(filename, R_OK) == -1)
    {
        printf("WARNING: reference file <%s> not found.\n", filename);
        return;
    }
    else
        printf("reading from <%s> ...\n", filename);

    std::string r = "r";
    FILE * f = fopen(filename, "r");

    assert(f);

    fread(&xmin_, sizeof(double), 1, f);
    fread(&xmax_, sizeof(double), 1, f);

    double M = 0;
    fread(&M, sizeof(double), 1 , f);
    M_ = int(M);
    nsources_ = M_*M_;

    posix_memalign((void **)&xsources_, 32, sizeof(double) * nsources_);
    posix_memalign((void **)&ysources_, 32, sizeof(double) * nsources_);
    posix_memalign((void **)&wsources_, 32, sizeof(double) * nsources_);

    fread(xsources_, sizeof(double), nsources_, f);
    fread(ysources_, sizeof(double), nsources_, f);
    fread(wsources_, sizeof(double), nsources_, f);

    posix_memalign((void **)&xmesh_, 32, sizeof(double) * M_*M_);
    posix_memalign((void **)&ymesh_, 32, sizeof(double) * M_*M_);

    posix_memalign((void **)&phi_, 32, sizeof(double) * M_*M_);

    posix_memalign((void **)&ux_, 32, sizeof(double) * M_*M_);
    posix_memalign((void **)&uy_, 32, sizeof(double) * M_*M_);

    posix_memalign((void **)&vorticity_, 32,sizeof(double) * M_*M_);
    posix_memalign((void **)&vorticity_temp_, 32, sizeof(double)* M_*M_);


}

double Vortexsimulation::getTime()
{
    return T_;
}

void Vortexsimulation::getMesh()
{

    //generate mesh from (xmin, ymin) to (xmax, ymax)
    for (int i=0; i<M_; ++i)
    {
        for (int j=0; j<M_; ++j)
        {
            xmesh_[i*M_+j] = xmin_ + j*dx_;
            ymesh_[i*M_+j] = xmin_ + i*dx_;
        }
    }
    celllist_.generateCelllist(cutoff_, dx_);
}

void Vortexsimulation::advance()
{
    timer t;
    t.start();
    poisson(xsources_, ysources_, wsources_, nsources_, phi_, M_, xmesh_, ymesh_, theta);
    t.stop();
    //std::cout << "poisson done in " << t.get_timing() << std::endl;
    t.start();
    velocityfield_openBC(phi_, M_, dx_, ux_, uy_);
    t.stop();
    //std::cout << "velocity done in " << t.get_timing() << std::endl;
    t.start();
    vorticityfield_by_curl_openBC(ux_, uy_, M_, dx_, vorticity_);
    t.stop();
    //std::cout << "vorticity done in " << t.get_timing() << std::endl;
    dt_ = getdt();
    t.start();
    particlexchangemethod_openBC(xmesh_, ymesh_, vorticity_, vorticity_temp_, M_, dx_, VISC, dt_, celllist_);
    t.stop();
    //std::cout << "PEM done in " << t.get_timing() << std::endl;
    t.start();
    convection(xsources_, ysources_, wsources_, ux_, uy_, xmesh_, ymesh_, vorticity_, M_, dt_);
    t.stop();
    T_ += dt_;
    //std::cout << "convection done in" << t.get_timing() << std::endl;
}

double Vortexsimulation::getdt()
{
    double max=0.0;
    for ( int i=0; i<M_*M_; ++i)
    {

        double tmp = ux_[i]*ux_[i] + uy_[i]*uy_[i];
        if(max<tmp && tmp<1000000.0 )
        {
            max=tmp;
        }
    }

    return LCFL/std::sqrt(max);

}

void Vortexsimulation::advanceandoutput(std::string f)
{

    timer t;
    t.start();
    poisson(xsources_, ysources_, wsources_, nsources_, phi_, M_, xmesh_, ymesh_, theta);
    t.stop();
    //std::cout << "poisson done in " << t.get_timing() << std::endl;
    t.start();
    velocityfield_openBC(phi_, M_, dx_, ux_, uy_);
    t.stop();
    //std::cout << "velocity done in " << t.get_timing() << std::endl;
    t.start();
    vorticityfield_by_curl_openBC(ux_, uy_, M_, dx_, vorticity_);
    t.stop();
    //std::cout << "vorticity done in " << t.get_timing() << std::endl;
    dt_ = getdt();

    std::string filenamecvs;
    ++outnumber_;

    std::string num= std::to_string(outnumber_);
    filenamecvs= "animationcvs/"+f+num;
    std::fstream cvsoutfile;
    cvsoutfile.open (filenamecvs, std::fstream::in | std::fstream::out | std::fstream::app);


    for (int i=0; i<nsources_; ++i)
    {

        cvsoutfile<< xmesh_[i] << " , "  << ymesh_[i] << " , "<< ux_[i] << " , "  << uy_[i] << " , "<<phi_[i] <<" , "<< vorticity_[i] <<" , "<< vorticity_[i]*1000.0 <<" , "<< wsources_[i] << std::endl;
    }

    cvsoutfile.close();

    t.start();
    particlexchangemethod_openBC(xmesh_, ymesh_, vorticity_, vorticity_temp_, M_, dx_, VISC, dt_, celllist_);
    t.stop();
    //std::cout << "PEM done in " << t.get_timing() << std::endl;
    t.start();
    convection(xsources_, ysources_, wsources_, ux_, uy_, xmesh_, ymesh_, vorticity_, M_, dt_);
    t.stop();
    T_ += dt_;


    //std::cout << "convection done in" << t.get_timing() << std::endl;

}

void Vortexsimulation::analyze()
{
    char filename[] = "./InitialConditions/vortex100_20t.dat";

    std::string r = "r";
    FILE * f = fopen(filename, "r");

    int nsources = 100*100;
    double *xmesh_analyt, *ymesh_analyt, *w_analyt, *w_error;
    posix_memalign((void **)&xmesh_analyt, 32, sizeof(double) * nsources);
    posix_memalign((void **)&ymesh_analyt, 32, sizeof(double) * nsources);
    posix_memalign((void **)&w_analyt, 32, sizeof(double) * nsources);
    posix_memalign((void **)&w_error, 32, sizeof(double) * nsources);

    fread(xmesh_analyt, sizeof(double), nsources, f);
    fread(ymesh_analyt, sizeof(double), nsources, f);
    fread(w_analyt, sizeof(double), nsources, f);

    std::string filenameh = "vortex100_error.dat";
    std::ofstream errorfile;
    errorfile.open(filenameh);
    errorfile << " X " << "\t"  << " Y "<< "\t" << " W "<< std::endl;

    for(int i = 0; i < nsources; ++i)
    {
        w_error[i] = std::abs(w_analyt[i]-vorticity_[i]);
        errorfile<< xmesh_[i] << "\t" << ymesh_[i] << "\t" << w_error[i] << std::endl;
    }

    errorfile.close();

}
