#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

//constants
#define pi (3.14159265)
#define c (30000.0 / Rs) // c in Rs/sec
#define omegaB (1.17 * 1.76 * pow(10.0, 19)) // e*Bo/(m*c)

//pulsar constants
#define P (1.382449) // period in sec
#define Rs (1.0) // R_star in 10 km
#define Omega (2.0 * pi / P) // 2 pi / P (same time z-component of Omega)
#define Rlc (c / Omega) // R_lc (light cylinder)
#define Rpc (Rs * sqrt(Rs / Rlc)) // R_pc (polar cap)

//model parameters
#define lambda (30000.0)
//#define f0 (0.5)
//#define A (0.1)
#define f0 (1.0)
#define A (5.0)
#define Rm (50.0 * Rs)

#define omega (2 * pi * 0.327 * pow(10.0, 9)) // 10^9 Hz = 1 GHz
#define gamma (100.0) // gamma factor

#define chi (48.0 * (pi / 180)) // chi in rad
#define beta (0.4 * (pi / 180)) // beta in rad

#define mode (1) // 0: X-mode, 1: O-mode

//-------------------------------------------------------------------------------
//functions
double r (double x, double y, double z); // r
double nr (double x, double y, double z, int a); // n_a
double CosPsim (double x, double y, double z, double mx, double my, double mz); // Cos(psi_m) = (m , n)
double rPerp (double x, double y, double z, double mx, double my, double mz); // r_perp
double g1 (double x, double y, double z, double mx, double my, double mz); // (f_0 R_0 / r_perp)^5
double g2 (double x, double y, double z, double mx, double my, double mz); // (R_0 / r_perp)^2
double gamma_a(double x, double y, double z, double mx, double my, double mz);  // g / x

double Ba (double x, double y, double z, double mx, double my, double mz, int a); // B_a
double B (double x, double y, double z, double mx, double my, double mz); // B
double ba (double x, double y, double z, double mx, double my, double mz, int a); // b

double g (double x, double y, double z, double mx, double my, double mz); // g(r_perp)
double theta (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz); // theta
double Theta (double x, double y, double z, double mx, double my, double mz); // Theta
double n (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz); // n

double Lambda (double x, double y, double z, double mx, double my, double mz); // Lambda

int Kdelta (int a, int i); // Kronecker Delta

//derivatives
// r_i derivatives
double dnadri (double x, double y, double z, int a, int i); // dn_a / dr_i
double mdndri (double x, double y, double z, double mx, double my, double mz, int i); // (m , dn / dr_i)
double dBadri (double x, double y, double z, double mx, double my, double mz, int a, int i); // dB_a / dr_i
double dbadri (double x, double y, double z, double mx, double my, double mz, int a, int i); // db_a / dr_i
double dthdri (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz, int i); // dth / dr_i
double drpdri (double x, double y, double z, double mx, double my, double mz, int i); // dr_perp / dr_i
double dg1dri (double x, double y, double z, double mx, double my, double mz, int i); // dg_1 / dr_i
double dg2dri (double x, double y, double z, double mx, double my, double mz, int i); // dg_2 / dr_i
double dgdri (double x, double y, double z, double mx, double my, double mz, int i); // dg / dr_i
double dThdri (double x, double y, double z, double mx, double my, double mz, int i); // dTh / dr_i

double dkndri (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz, int i); // d(k/n)/dr_i

// k_i derivatives
double dkadki (double kx, double ky, double kz, int a, int i); // d(k_a/k)/dk_i
double dkndki (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz, int i); // d(k/n)/dk_i

void showVector(const char *msg, double x, double y, double z) {
  cout << msg << "x: " << x << " | "
       << msg << "y: " << y << " | "
       << msg << "z: " << z << "\n";
}

int main()
{
    ofstream par;
 //   par.open ("data.dat");

    // geometry of observing surface
    double ox = sin(chi + beta);
    double oy = 0;
    double oz = cos(chi + beta);
    double p1x = cos(chi + beta);
    double p1y = 0;
    double p1z = -sin(chi + beta);
    double p2x = 0;
    double p2y = 1;
    double p2z = 0;

    // intensity
    double I;
    double sumInt;

    // initial phase
    double fi0;

    // a-parameters
    double a1, a2;
    double a1min = -10;
    double a1max = 10; // this is the perpendicular to movement coordinate
    double a2min = -10;
    double a2max = 10; // this is along the movement
    double astep = 0.25;

    // fi-parameters
    double fi;
    double fimax = 15;
    double fistep = 0.5;

    double Rmax = 1000 * Rs;
    double Tmax = Rmax / c;
    double rx, ry, rz;
    double mx, my, mz = cos(chi);

    double kx, ky, kz;

    double t, dt = 0.00005; // normally 0.00005

    // counter
    double counter = 0;
    for(fi = -fimax; fi <= fimax; fi += fistep)
        for(a1 = a1min; a1 <= a1max; a1 += astep)
            for(a2 = a2min; a2 <= a2max; a2 += astep)
                counter ++;
    double cmax = counter;
    counter = 0;

    int run_count = 0; //andrey
    for(fi = -fimax; fi <= fimax; fi += fistep) // cycle for fi
    {
/*  andrey */
        std::string iter_str = std::to_string(run_count);
        int nzeros = 0;
        if (run_count < 10) {
          nzeros = 2;
        } else if (run_count < 100) {
          nzeros = 1;
        }
        run_count++;
        iter_str = std::string(nzeros, '0').append(iter_str);
        string s_file = "tmp/output/xmode_";
        s_file = s_file + iter_str;//to_string(run_count);
        s_file = s_file + "_phi_";
        s_file = s_file + to_string(fi);
        s_file = s_file + ".dat";
        par.open (s_file); // andrey
/*  andrey */

        sumInt = 0;
        		// cout << fi << endl << endl; // counter
        fi0 = fi * pi / 180;
        cout << (int)(counter*100/cmax) << " %" << endl; // counter
        int temppppp = 0;
        for(a1 = a1min; a1 <= a1max; a1 += astep) // cycle for a1
        {
            cout << (counter*100/cmax) << " %" << endl; // counter
            for(a2 = a2min; a2 <= a2max; a2 += astep) //
            {
                counter ++;
                I = 0;
                rx = Rmax*ox + a1*p1x + a2*p2x;
                ry = Rmax*oy + a1*p1y + a2*p2y;
                rz = Rmax*oz + a1*p1z + a2*p2z;
                kx = -ox;
                ky = -oy;
                kz = -oz;
                for(t = Tmax; t >= -Tmax / 10; t -= dt)
                {
                    if(r(rx,ry,rz) < 1.2*Rs)
                        break;

                    mx = sin(chi)*cos(fi0 + Omega*t);
                    my = sin(chi)*sin(fi0 + Omega*t);

                    if(theta(rx,ry,rz,mx,my,mz,kx,ky,kz) > pi/2.)
                        break;

                    rx += c*dkndki(rx,ry,rz,mx,my,mz,kx,ky,kz,1)*dt;
                    ry += c*dkndki(rx,ry,rz,mx,my,mz,kx,ky,kz,2)*dt;
                    rz += c*dkndki(rx,ry,rz,mx,my,mz,kx,ky,kz,3)*dt;

                    kx += -c*dkndri(rx,ry,rz,mx,my,mz,kx,ky,kz,1)*dt;
                    ky += -c*dkndri(rx,ry,rz,mx,my,mz,kx,ky,kz,2)*dt;
                    kz += -c*dkndri(rx,ry,rz,mx,my,mz,kx,ky,kz,3)*dt;

                    mx = sin(chi)*cos(fi0 + Omega*(t-dt));
                    my = sin(chi)*sin(fi0 + Omega*(t-dt));

                    I += g(rx,ry,rz,mx,my,mz)*exp(-A*pow((r(rx,ry,rz) - Rm),2)/(Rm*Rm))*exp(-pow(theta(rx,ry,rz,mx,my,mz,kx,ky,kz)*gamma_a(rx,ry,rz,mx,my,mz),2));
                }
                par << a1 << "," << a2 << "," << I << endl;
            }
        }
        par << endl;
        par.close(); //andrey
    }
//    par.close();
    return 0;
}

//functions
double r (double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}
double nr (double x, double y, double z, int a)
{
    if (a == 1)
        return x/r(x,y,z);
    if (a == 2)
        return y/r(x,y,z);
    else
        return z/r(x,y,z);
}
double CosPsim (double x, double y, double z, double mx, double my, double mz)
{
    double nx = nr(x,y,z,1);
    double ny = nr(x,y,z,2);
    double nz = nr(x,y,z,3);
    return nx*mx + ny*my + nz*mz;
}
double rPerp (double x, double y, double z, double mx, double my, double mz)
{
    double R = r (x,y,z);
    double cospsim = CosPsim (x,y,z,mx,my,mz);
    if (1 - pow(cospsim,2) < 0)
        cout << "ERROR!" << endl;
    return sqrt(1 - pow(cospsim,2))*sqrt(Rs/R)*Rs;
}
double g1 (double x, double y, double z, double mx, double my, double mz)
{
    double rperp = rPerp (x,y,z,mx,my,mz);
    return pow(sqrt(f0)*Rpc/rperp,2);
}
double g2 (double x, double y, double z, double mx, double my, double mz)
{
    double rperp = rPerp (x,y,z,mx,my,mz);
    return pow(Rpc/rperp,2);
}

double Ba (double x, double y, double z, double mx, double my, double mz, int a)
{
    double nx = nr (x,y,z,1);
    double ny = nr (x,y,z,2);
    double nz = nr (x,y,z,3);
    double R3 = pow(r(x,y,z),3);
    double nn[3] = {nx,ny,nz};
    double mm[3] = {mx,my,mz};
    return 3*CosPsim(x,y,z,mx,my,mz)*nn[a - 1]/R3 - mm[a - 1]/R3;
}
double B (double x, double y, double z, double mx, double my, double mz)
{
    double Bx = Ba (x,y,z,mx,my,mz,1);
    double By = Ba (x,y,z,mx,my,mz,2);
    double Bz = Ba (x,y,z,mx,my,mz,3);
    return sqrt(Bx*Bx + By*By + Bz*Bz);
}
double ba (double x, double y, double z, double mx, double my, double mz, int a)
{
    double Bx = Ba (x,y,z,mx,my,mz,1);
    double By = Ba (x,y,z,mx,my,mz,2);
    double Bz = Ba (x,y,z,mx,my,mz,3);
    double BB = B (x,y,z,mx,my,mz);
    double Bb[3] = {Bx,By,Bz};
    return Bb[a - 1]/BB;
}

double g (double x, double y, double z, double mx, double my, double mz)
{
    double rperp = rPerp(x, y, z, mx, my, mz);
    double cospsim = CosPsim(x,y,z,mx,my,mz);
    double sinpsim = sqrt(1.0 - cospsim * cospsim);
    double rperp_over_R0 = sinpsim * sqrt(Rlc / r(x, y, z));
//    return exp(-pow(rperp_over_R0, 2)) / (1.0 + pow(f0 / rperp_over_R0, 5));
//    return exp(-pow(2*rperp_over_R0, 2)) / (1.0 + pow(f0 / rperp_over_R0, 2)); //andrey
//    return 1.25 * sqrt(sqrt(rperp_over_R0*rperp_over_R0) * 2.5) * exp(-pow(rperp_over_R0*2.5,2));//andrey
    double K = 14.8;

//    return K * pow(rperp_over_R0, 3) * exp(-10 * pow(rperp_over_R0,4));//andrey
    //return 100.0/1.284 * pow(rperp_over_R0, 3) * exp(-10 * pow(rperp_over_R0,2));//andrey
//    return 130.0 * pow(rperp_over_R0,3) * exp(-10.0 * pow(rperp_over_R0,5)) * pow(1.0 - rperp_over_R0,2);

   return (37.5 * pow(rperp_over_R0,3) * exp(-24.0 * pow(rperp_over_R0,6) ));

}

double gamma_a (double x, double y, double z, double mx, double my, double mz)
{
    double rperp = rPerp(x,y,z,mx,my,mz);
    double cospsim = CosPsim(x,y,z,mx,my,mz);
    double sinpsim = sqrt(1.0 - cospsim * cospsim);
    double rperp_over_R0 = sinpsim * sqrt(Rlc / r(x,y,z));

    if ((rperp_over_R0) < 1e-10){rperp_over_R0 = 1e-10;}
    return gamma / rperp_over_R0;
}

double theta (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz)
{
    double bx = ba(x,y,z,mx,my,mz,1);
    double by = ba(x,y,z,mx,my,mz,2);
    double bz = ba(x,y,z,mx,my,mz,3);
    double th = acos(-(bx*kx + by*ky + bz*kz)/sqrt(kx*kx + ky*ky + kz*kz));
    return th;
}
double Theta (double x, double y, double z, double mx, double my, double mz)
{
    double G = g (x,y,z,mx,my,mz);
    double R = r (x,y,z);
    double bz = ba(x,y,z,mx,my,mz,3);
    double Lamb = Lambda(x,y,z,mx,my,mz);
    return fabs(Lamb*lambda*G/pow(R,3));
}
double n (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz)
{
    double th = theta (x,y,z,mx,my,mz,kx,ky,kz);
    double Th = Theta (x,y,z,mx,my,mz);
    return 1 + th*th/4 - sqrt(pow(th,4)/16 + Th);
}

double Lambda (double x, double y, double z, double mx, double my, double mz)
{
    double R = r(x,y,z);
    double bz = ba(x,y,z,mx,my,mz,3);
    return 2*omegaB*Omega*bz/(omega*omega*gamma_a(x,y,z,mx,my,mz)*gamma_a(x,y,z,mx,my,mz)*gamma_a(x,y,z,mx,my,mz));
}

int Kdelta (int a, int i)
{
    if (a == i)
        return 1;
    else
        return 0;
}

//derivatives
// r_i derivatives
double dnadri (double x, double y, double z, int a, int i)
{
    double R = r (x,y,z);
    double rvec[3] = {x,y,z};
    return (Kdelta(a,i) - nr(x,y,z,a)*nr(x,y,z,i))/R;
}
double mdndri (double x, double y, double z, double mx, double my, double mz, int i)
{
    return	mx*dnadri(x,y,z,1,i) + my*dnadri(x,y,z,2,i) + mz*dnadri(x,y,z,3,i);
}
double dBadri (double x, double y, double z, double mx, double my, double mz, int a, int i)
{
    double R = r(x,y,z);
    double Mn = CosPsim (x,y,z,mx,my,mz);
    double mDn = mdndri (x,y,z,mx,my,mz,i);
    double nvec[3] = {x/R,y/R,z/R};
    double Mvec[3] = {mx,my,mz};
    double rvec[3] = {x,y,z};

    return - 3*rvec[i - 1]*(3*Mn*nvec[a - 1] - Mvec[a - 1])/pow(R,5) + (3/pow(R,3))*dnadri(x,y,z,a,i)*Mn + (3/pow(R,3))*nvec[a - 1]*mDn;
}
double dbadri (double x, double y, double z, double mx, double my, double mz, int a, int i)
{
    double Bvec[3] = {Ba(x,y,z,mx,my,mz,1),Ba(x,y,z,mx,my,mz,2),Ba(x,y,z,mx,my,mz,3)};
    double dBvec[3] = {dBadri(x,y,z,mx,my,mz,1,i), dBadri(x,y,z,mx,my,mz,2,i), dBadri(x,y,z,mx,my,mz,3,i)};
    double BB = B (x,y,z,mx,my,mz);
    return (1/BB)*dBvec[a - 1] - (Bvec[a - 1]/(BB*BB*BB))*(Bvec[0]*dBvec[0] + Bvec[1]*dBvec[1] + Bvec[2]*dBvec[2]);
}

double drpdri (double x, double y, double z, double mx, double my, double mz, int i)
{
    double cos = fabs(CosPsim (x,y,z,mx,my,mz)), sin = sqrt(1 - cos*cos);
    if (cos > 1)
        cout << "ERROR!" << endl;
    double R = r (x,y,z), rvec[3] = {x,y,z};
    double mdn = mdndri (x,y,z,mx,my,mz,i);
    return - Rs*sqrt(Rs)*(rvec[i - 1]*sin/(2*R*R) + cos*mdn/sin)/sqrt(R);
}
double dg1dri (double x, double y, double z, double mx, double my, double mz, int i)
{
    double drp = drpdri (x,y,z,mx,my,mz,i);
    double g01 = g1(x,y,z,mx,my,mz);
    double rp = rPerp(x,y,z,mx,my,mz);
    return -5*g01*drp/rp;
}
double dg2dri (double x, double y, double z, double mx, double my, double mz, int i)
{
    double g02 = g2(x,y,z,mx,my,mz);
    double rp = rPerp(x,y,z,mx,my,mz);
    double drp = drpdri(x,y,z,mx,my,mz,i);
    return -2*g02*drp/rp;
}

double andrey_dgdri (double x, double y, double z, double mx, double my, double mz, int i)
{
//    double Dg1 = dg1dri (x,y,z,mx,my,mz,i);
//    double Dg2 = dg2dri (x,y,z,mx,my,mz,i);
//    double G = g (x,y,z,mx,my,mz);
//    double g01 = g1 (x,y,z,mx,my,mz);
//    double g02 = g2 (x,y,z,mx,my,mz);
//    double r_i = 1;
//    double A = 1;
//    double B = 1;
    return 1;//G*((1 + g01)*Dg2 - g02*g02*Dg1)/(g02*g02*(1 + g01));
}


double dgdri (double x, double y, double z, double mx, double my, double mz, int i)
{
    double Dg1 = dg1dri (x,y,z,mx,my,mz,i);
    double Dg2 = dg2dri (x,y,z,mx,my,mz,i);
    double G = g (x,y,z,mx,my,mz);
    double g01 = g1 (x,y,z,mx,my,mz);
    double g02 = g2 (x,y,z,mx,my,mz);
    return G*((1 + g01)*Dg2 - g02*g02*Dg1)/(g02*g02*(1 + g01));
//    return andrey_dgdri (x,y,z,mx,my,mz);
}
double dThdri (double x, double y, double z, double mx, double my, double mz, int i)
{
    double rvec[3] = {x,y,z};
    double G = g (x,y,z,mx,my,mz);
    double R = r (x,y,z);
    double DgDri = dgdri(x,y,z,mx,my,mz,i);
    double Lamb = Lambda (x,y,z,mx,my,mz);
    return Lamb*lambda*pow(Rs,3)*(1/pow(R,3))*(dgdri(x,y,z,mx,my,mz,i) - (3*rvec[i - 1]/(R*R))*G);
}

double dkndri (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz, int i)
{
    if (mode == 0)
        return 0; //
    double k = sqrt (kx*kx + ky*ky + kz*kz);
    double N = n (x,y,z,mx,my,mz,kx,ky,kz);
    double bk = kx*dbadri(x,y,z,mx,my,mz,1,i) + ky*dbadri(x,y,z,mx,my,mz,2,i) + kz*dbadri(x,y,z,mx,my,mz,3,i);
    double th = theta (x,y,z,mx,my,mz,kx,ky,kz);
    double Th = Theta (x,y,z,mx,my,mz);
    double DThdri = dThdri (x,y,z,mx,my,mz,i);
    if (pow(th,4)/16+Th < 0)
        cout << "ERROR!" << endl;
    return k*(bk - (pow(th,2)*bk/4 - DThdri)/sqrt(pow(th,4)/16+Th))/(2*N*N);
}

// k_i derivatives
double dkadki (double kx, double ky, double kz, int a, int i)
{
    double K = sqrt(kx*kx + ky*ky + kz*kz);
    double kvec[3] = {kx,ky,kz};
    return Kdelta(a,i)/K - kvec[a - 1]*kvec[i - 1]/pow(K,3);
}
double dkndki (double x, double y, double z, double mx, double my, double mz, double kx, double ky, double kz, int i)
{
    double kvec[3] = {kx,ky,kz};
    if (mode == 0)
        return kvec[i - 1]; //
    double dkvec[3] = {dkadki(kx,ky,kz,1,i),dkadki(kx,ky,kz,2,i),dkadki(kx,ky,kz,3,i)};
    double N = n (x,y,z,mx,my,mz,kx,ky,kz);
    double th = theta (x,y,z,mx,my,mz,kx,ky,kz);
    double Th = Theta (x,y,z,mx,my,mz);
    double bvec[3] = {ba (x,y,z,mx,my,mz,1), ba (x,y,z,mx,my,mz,2), ba (x,y,z,mx,my,mz,3)};
    double K = sqrt(kx*kx + ky*ky + kz*kz);
    return kvec[i - 1]/(N*K) + (K/(2*N*N))*(1 - (th*th/4)/sqrt(pow(th,4)/16 + Th))*(bvec[0]*dkvec[0] + bvec[1]*dkvec[1] + bvec[2]*dkvec[2]);
}
