/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

struct UserData
{
    double sawL;
    double g;
    double theta;
    double Tf;
};

void AddSawPlate(DEM::Domain & dom, int Tag,  Vec3_t & X, double Lx, double Ly, size_t Ntooth, double depth, double rho, double R)
{
     Array<Vec3_t>        V(4*Ntooth+2);
     Array<Array <int> >  E(2*Ntooth+1);
     Array<Array <int> >  F(2*Ntooth);
     double step = Lx/Ntooth;
     for (size_t i=0;i<Ntooth+1;i++)
     {
         V[i             ] = Vec3_t(i*step,-0.5*Ly,0.0);
         V[i + Ntooth + 1] = Vec3_t(i*step, 0.5*Ly,0.0);
         E[i].Push(i);
         E[i].Push(i+Ntooth+1);
     }
     //std::cout << "1" << std::endl;
     for (size_t i=0;i<Ntooth;i++)
     {
         V[i + 2*Ntooth + 2] = Vec3_t(i*step+0.5*step,-0.5*Ly,depth);
         V[i + 3*Ntooth + 2] = Vec3_t(i*step+0.5*step, 0.5*Ly,depth);
         E[i + Ntooth + 1].Push(i + 2*Ntooth + 2);
         E[i + Ntooth + 1].Push(i + 3*Ntooth + 2);
         F[2*i  ].Push(i               );
         F[2*i  ].Push(i +   Ntooth + 1);
         F[2*i  ].Push(i + 3*Ntooth + 2);
         F[2*i  ].Push(i + 2*Ntooth + 2);
         F[2*i+1].Push(i + 1           );
         F[2*i+1].Push(i +   Ntooth + 2);
         F[2*i+1].Push(i + 3*Ntooth + 2);
         F[2*i+1].Push(i + 2*Ntooth + 2);
     }
     //std::cout << "2" << std::endl;
     dom.Particles.Push(new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
     //std::cout << "3" << std::endl;
     dom.Particles[dom.Particles.Size()-1]->Q          = 1.0,0.0,0.0,0.0;
     dom.Particles[dom.Particles.Size()-1]->Props.V    = Lx*Ly*R;
     dom.Particles[dom.Particles.Size()-1]->Props.m    = rho*Lx*Ly*R;
     dom.Particles[dom.Particles.Size()-1]->I          = 1.0,1.0,1.0;
     dom.Particles[dom.Particles.Size()-1]->I         *= dom.Particles[dom.Particles.Size()-1]->Props.m;
     dom.Particles[dom.Particles.Size()-1]->x          = Vec3_t(0.5*Lx,0.0,depth);
     dom.Particles[dom.Particles.Size()-1]->Ekin       = 0.0;
     dom.Particles[dom.Particles.Size()-1]->Erot       = 0.0;
     dom.Particles[dom.Particles.Size()-1]->Dmax       = sqrt(Lx*Lx+Ly*Ly+depth*depth)+R;
     dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
     dom.Particles[dom.Particles.Size()-1]->Index      = dom.Particles.Size()-1;

     dom.Particles[dom.Particles.Size()-1]->Position(X);
}

void Setup (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    
    double scale = std::min(100.0*dom.Time/dat.Tf,1.0);
    for (size_t np=0;np<dom.Particles.Size();np++)
    {
        dom.Particles[np]->Ff = dom.Particles[np]->Props.m*Vec3_t(dat.g*sin(scale*dat.theta*M_PI/180.0),0.0,-dat.g*cos(scale*dat.theta*M_PI/180.0));
    }
}

void Report (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    
    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    size_t  RenderVideo;// Decide is video should be render
    bool   Cohesion;    // Decide if coheison is going to be simulated
    double fraction;    // Fraction of particles to be generated
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Beta;        // Rolling stiffness coefficient (only for spheres)
    double Eta;         // Plastic moment coefficient (only for spheres, 0 if rolling resistance is not used)
    double Eps;         // Threshold for breking bonds
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    double g;           // gravity
    double theta;       // slope
    double T0;          // Time span for the compression
    double Tf;          // Final time for the test
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Cohesion;     infile.ignore(200,'\n');
        infile >> fraction;     infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Beta;         infile.ignore(200,'\n');
        infile >> Eta;          infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> g;            infile.ignore(200,'\n');
        infile >> theta;        infile.ignore(200,'\n');
        infile >> T0;           infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
    }

    // domain and User data
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha=verlet;
    dom.Dilate = true;

    bool load = false;
    // particle
    if      (ptype=="sphere")  dom.GenSpheres  (-1, Lx, nx, rho, "HCP", seed, fraction);
    else if (ptype=="sphereboxhcp") 
    {
        Vec3_t Xmin(-0.5*Lx,-0.5*Ly,-0.5*Lz);
        Vec3_t Xmax = -Xmin;
        dom.GenSpheresBox (-1, Xmin, Xmax, R, rho, "HCP",    seed, fraction, Eps);
    }
    else if (ptype=="voronoi")
    {
        if (ny==1) dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, Cohesion, bVec3_t(true,false,true), seed, fraction, Vec3_t(0.0,1.0,0.0));
        else       dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, Cohesion, bVec3_t(true,true ,true), seed, fraction, Vec3_t(0.0,0.0,0.0));
    }
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh,/*R*/R,/*rho*/rho,Cohesion,false);
    }
    else if (ptype=="rice") dom.GenRice(-1,Lx,nx,R,rho,seed,fraction);
    else 
    {
        dom.Load(ptype.CStr());
        load = true;
    }

    if (!load)
    {    
        Vec3_t Xmin,Xmax;
        dom.BoundingBox(Xmin,Xmax);

        Vec3_t X0(0.0,Xmax(1)+2*R,0.0);
        Vec3_t X1(0.0,Xmin(1)-2*R,0.0);
        dom.AddPlane(-2,X0,R,1.3*Lx,1.3*Lz,3.0,0.5*M_PI,&OrthoSys::e0);
        dom.AddPlane(-3,X1,R,1.3*Lx,1.3*Lz,3.0,0.5*M_PI,&OrthoSys::e0);

        X0 = Vec3_t(0.0,0.0,Xmin(2)-2*R);
        X1 = Vec3_t(0.0,0.0,Xmax(2)+2*R);
        //dom.AddPlane(-4,X0,R,1.3*Lx,1.3*Ly,3.0,0.0,&OrthoSys::e0);
        //dom.AddPlane(-5,X1,R,1.3*Lx,1.3*Ly,3.0,0.0,&OrthoSys::e0);
        AddSawPlate(dom,-4,X0,Lx*2.0,2.0*Ly,8,Lx/4.0,3.0,R);
        //AddSawPlate(dom,-5,X1,Lx*2.0,2.0*Ly,8,Lx/4.0,3.0,R);
        //Quaternion_t q;
        //NormalizeRotation (M_PI,OrthoSys::e1,q);
        //dom.GetParticle(-4)->Rotate(q,dom.GetParticle(-5)->x);
    }



    dom.Xmax =  0.5*Lx;
    dom.Xmin = -0.5*Lx;

    Dict B1;
    B1.Set(-1,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt ,Mu ,Beta,Eta);
    B1.Set(-2,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,0.0,0.0,Beta,Eta);
    B1.Set(-3,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,0.0,0.0,Beta,Eta);
    B1.Set(-4,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,0.0,Mu ,Beta,Eta);
    dom.SetProps(B1);

    dom.GetParticle(-2)->FixVeloc();
    dom.GetParticle(-3)->FixVeloc();
    dom.GetParticle(-4)->FixVeloc();

    String fkey_a  (filekey+"_a");
    String fkey_b  (filekey+"_b");
    String fkeybf_a(filekey+"bf_a");
    String fkeybf_b(filekey+"bf_b");

    for (size_t np=0;np<dom.Particles.Size();np++)
    {
        dom.Particles[np]->Ff = dom.Particles[np]->Props.m*Vec3_t(0.0,0.0,-g);
    }

    dt = dt*dom.CriticalDt();
    dom.Solve (/*tf*/T0, /*dt*/dt, /*dtOut*/dtOut, NULL, NULL, fkey_a.CStr(),RenderVideo,Nproc);
    dom.Save(fkey_a.CStr());
    dom.WriteXDMF(fkey_a.CStr());
    dom.WriteBF(fkeybf_a.CStr());

    dat.sawL = Lx/4.0;
    dat.g    = g;
    dat.theta= theta;
    dat.Tf   = Tf;

    dom.Time = 0.0;
    dom.Solve (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_b.CStr(),RenderVideo,Nproc);
    dom.Save(fkey_b.CStr());
    dom.WriteXDMF(fkey_b.CStr());
    dom.WriteBF(fkeybf_b.CStr());

    return 0;
}
MECHSYS_CATCH
