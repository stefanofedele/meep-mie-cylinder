#include <meep.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "materials.h"
#include "functions.h"
#include "parameters.h"

using namespace std;
using namespace meep;

//simulation parameters
parameters m;
const double lx = m.v["lx"], ly = m.v["ly"];
const double cyl_rad = m.v["cyl_rad"];
const double cyl_eps = m.v["cyl_eps"];
const double cyl_conduc = m.v["cyl_conduc"];
const int num_freq = m.v["num_freq"];
const double freq = m.v["freq"], fwidth = m.v["fwidth"];
const double f_min = m.v["f_min"], f_max = m.v["f_max"];
const double amicron = m.v["amicron"];
const double resolution = m.v["resolution"];

//defines a single material location function
bool where(const vec &p) {
    if (pow((p.x()-lx/2),2) + pow(p.y()-ly/2,2) < pow(cyl_rad,2))
        return true;
    return false;
}

double vac(const vec &p) {return 1.0;}

double eps(const vec &p) {
    if (where(p)) return cyl_eps;
    else return 1.0;
}

double conduc(const vec &p) {
    return cyl_conduc*where(p);
}

int main(int argc, char **argv) {
    initialize mpi(argc, argv);
    grid_volume v = vol2d(lx,ly,resolution);
    structure sf(v, vac, pml(1));
    structure sg(v, eps, pml(1));
    //sg.set_conductivity(Ex,conduc);
    //material silver = Au(amicron*1e-6);
    //add_material(silver,where,sg); 

    fields f(&sf),g(&sg);
    f.use_real_fields();  //???
    g.use_real_fields();

    h5file eps_file("./output/eps.h5",h5file::WRITE,false);
    g.output_hdf5(Dielectric, v.surroundings(),&eps_file);

//Add the sources
    //gaussian_src_time src(freq, fwidth);
    //gaussian_src_time src2(freq*1.5,fwidth);
    //continuous_src_time src(2.0375,10);
    //volume src_vol(vec(lx/2.0,1),vec(lx/2.0,1));

    //volume src_vol(vec(1,1),vec(lx-1,1));
    volume src_vol(vec(lx/2,ly/2-cyl_rad-1),vec(lx/2,ly/2-cyl_rad-1));
    //gaussian_src_time src(freq, fwidth);
    continuous_src_time src(freq,.1);
    f.add_volume_source(Ex, src, src_vol);
    g.add_volume_source(Ex, src, src_vol);
    //f.add_volume_source(Ex, src2, src_vol);
    //g.add_volume_source(Ex, src2, src_vol);

//Add the flux regions    
    double buff = 1.0;
    double diag = cyl_rad*buff;
    vec center(lx/2, ly/2);
    vec lb(center + vec(-diag,-diag)), rb(center + vec(diag,-diag));
    vec lt(center + vec(-diag,diag)), rt(center + vec(diag,diag));
    vec li(center + vec(-diag,0)), ri(center + vec(diag,0));

    //wrap my volumes/flux regions into a class
    volume left(lb, lt);
    volume right(rb, rt);
    volume top(lt, rt);
    volume bottom(lb, rb);
    volume inc(li,ri);

    dft_flux fl = f.add_dft_flux_plane(left,f_min,f_max,num_freq);
    dft_flux fr = f.add_dft_flux_plane(right,f_min,f_max,num_freq);
    dft_flux ft = f.add_dft_flux_plane(top,f_min,f_max,num_freq);
    dft_flux fb = f.add_dft_flux_plane(bottom,f_min,f_max,num_freq);
    dft_flux finc = f.add_dft_flux_plane(inc,f_min,f_max,num_freq);

    dft_flux fl_g = g.add_dft_flux_plane(left,f_min,f_max,num_freq);
    dft_flux fr_g = g.add_dft_flux_plane(right,f_min,f_max,num_freq);
    dft_flux ft_g = g.add_dft_flux_plane(top,f_min,f_max,num_freq);
    dft_flux fb_g = g.add_dft_flux_plane(bottom,f_min,f_max,num_freq);

//creates own output h5 file. WRITE defined in h5file class. false = no parallel
    h5file field_file("./output/data.h5",h5file::WRITE,false);
    h5file field_file_norm("./output/data_norm.h5",h5file::WRITE,false);

//do normalization simulation
    int count = 0;
    while (f.time() < 80) {
        f.step();
        if (count == 50) {
            f.output_hdf5(Ex, v.surroundings(),&field_file_norm,true);
            count = 0;
        }
        count ++;
    }

//subtract off the incident for scattering simulation
    fl_g -= fl;
    fr_g -= fr;
    ft_g -= ft;
    fb_g -= fb;

//do scattering simulation
    count = 0;
    while (g.time() < 80) {
        g.step();
        if (count == 50) {
            g.output_hdf5(Ex, v.surroundings(),&field_file,true);
            count = 0;
        }
        count ++;
    }

    //Note:Maybe there's a list of dft_fluxes in fields that allows this to be a loop
//compute and output scattering fluxes
    double* flux_l = fl_g.flux();
    double* flux_r = fr_g.flux();
    double* flux_t = ft_g.flux();
    double* flux_b = fb_g.flux();
    double* flux_inc = finc.flux();

    output_data(flux_l, num_freq, "output/flux_l.dat");
    output_data(flux_r, num_freq, "output/flux_r.dat");
    output_data(flux_t, num_freq, "output/flux_t.dat");
    output_data(flux_b, num_freq, "output/flux_b.dat");
    output_data(flux_inc, num_freq, "output/flux_inc.dat");


    return 0;
}
