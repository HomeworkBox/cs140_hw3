/*
 Assignment 3
 Team Member 1 :
 Team Member 2 :
 */

#include "nBody.h"

#define G 6.674e-11

void readnbody(double* s, double* v, double* m, int n) {
    int myrank;
    int nprocs;
    int i,j;
    int bpp;
    double *tmpS, *tmpV, *tmpM;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Status status;
    bpp = n/nprocs;
    
    tmpS = (double *)malloc(sizeof(double) * bpp * 3);
    tmpV = (double *)malloc(sizeof(double) * bpp * 3);
    tmpM = (double *)malloc(sizeof(double) * bpp * 3);
    
    
    if (myrank == 0) {
        for (i = 0; i < bpp; i++) {
            double x, y, z, vx, vy, vz, mm;
            
            int result = scanf(INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mm);
            if (result != 7) {
                fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
                exit(0);
            }
            s[i * 3] = x;
            s[i * 3 + 1] = y;
            s[i * 3 + 2] = z;
            v[i * 3] = vx;
            v[i * 3 + 1] = vy;
            v[i * 3 + 2] = vz;
            m[i] = mm;
        }
        for (j = 1; j < nprocs; j++) {
            for (i = 0; i < n / nprocs; i++) {
                double x, y, z, vx, vy, vz, mm;
                
                int result = scanf(INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mm);
                if (result != 7) {
                    fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
                    exit(0);
                }
                tmpS[i * 3] = x;
                tmpS[i * 3 + 1] = y;
                tmpS[i * 3 + 2] = z;
                tmpV[i * 3] = vx;
                tmpV[i * 3 + 1] = vy;
                tmpV[i * 3 + 2] = vz;
                tmpM[i] = mm;
            }
            MPI_Send(tmpS, bpp * 3, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            MPI_Send(tmpV, bpp * 3, MPI_DOUBLE, j, 1, MPI_COMM_WORLD);
            MPI_Send(tmpM, bpp, MPI_DOUBLE, j, 2, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(s, bpp * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&status);
        MPI_Recv(v, bpp * 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,&status);
        MPI_Recv(m, bpp, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,&status);
    }
    
    free(tmpS);
    free(tmpV);
    free(tmpM);
}

void gennbody(double* s, double* v, double* m, int n) {
    printf("Generate nBody initial condition here.\n");
    if (n == 0) {
        m[0] = 1.99e30;
        m[1] = 3.30e23;
        m[2] = 4.87e24;
        m[3] = 5.97e24;
        m[4] = 6.42e23;
        m[5] = 1.90e27;
        m[6] = 5.69e26;
        m[7] = 1.02e26;
        m[8] = 8.68e25;
        m[9] = 1.31e22;
        double r[10] = {0, 5.8e10, 1.08e11, 1.50e11, 2.28e11, 7.78e11, 1.43e12, 2.87e12, 4.50e12, 5.91e12};
        int h;

        double vcopy[30] = {0.0000e+00,   0.0000e+00,   0.0000e+00,
            -4.0326e+04,   2.5893e+04   ,0.0000e+00,
            -2.9552e+04 ,  1.8975e+04   ,0.0000e+00,
            -2.5076e+04  , 1.6101e+04,   0.0000e+00,
            -2.0339e+04  , 1.3060e+04 ,  0.0000e+00,
            -1.1011e+04  , 7.0698e+03  , 0.0000e+00,
            -8.1214e+03  , 5.2147e+03,   0.0000e+00,
            -5.7327e+03  , 3.6809e+03 ,  0.0000e+00,
            -4.5782e+03  , 2.9396e+03  , 0.0000e+00,
            -3.9949e+03   ,2.5651e+03   ,0.0000e+00};
        for (h = 0; h < 10; h++){
            s[h*3] = r[h]*cos(1);
            s[h*3+1] = r[h]*sin(1);
            s[h*3+2] = 0;
            v[h*3] = vcopy[h*3];
            v[h*3+1] = vcopy[h*3+1];
            v[h*3+2] = vcopy[h*3+2];
            
        }
        
    }
    else if (n > 0) {
        double dist;
        double theta;
        int i;
        int p;
        int bpp;
        MPI_Comm_size(MPI_COMM_WORLD, &p);
        bpp = n/p;
        srand(time(NULL));
        
        for(i = 0; i < bpp; i++) {
            m[i] = ((double)rand() / (double)RAND_MAX) * 1e30;
            dist = ((double)rand() / (double)RAND_MAX) * 5e13;
            theta = ((double)rand() / (double)RAND_MAX) * (2 * 3.1415926);
            s[i * 3] = dist*cos(theta);
            s[i * 3 + 1] = dist*sin(theta);
            s[i * 3 + 2] = (((double)rand() / (double)RAND_MAX) - 0.5) * 1e11;
        }
    }
}

void nbody(double* s, double* v, double* m, int n, int iter, int timestep) {
    int myrank;
    int nprocs;
    int i, it;
    double* ss, *vv, *tmpPtr;
    MPI_Status status;
    int TAG3 = 3, TAG4 = 4, TAG5 = 5, TAG6 = 6; // tags for MPI_Send and MPI_RECV
    if (n == 0) {
        n = 10;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int bpp = n/nprocs;
    // if only 1 processor, data sent to itself, so TAG5 and TAG6 should be same to TAG3 and TAG4
    if (nprocs == 1) {
        TAG5 = 3;
        TAG6 = 4;
    }

    
    timestep = timestep * 604800; // Step is a week in seconds
    
    ss = (double *)malloc(sizeof(double) * bpp * 3); // MGR1: ss, an array of body position coordinates, size: bpp*3
    tmpPtr = ss; // pointer helps reduce number of mutiplications
    for (i = 0; i < bpp; i++) {
        *(tmpPtr++) = s[i * 3];
        *(tmpPtr++) = s[i * 3 + 1];
        *(tmpPtr++) = s[i * 3 + 2];
    }
    
    vv = (double *)malloc(sizeof(double) * bpp * 3);
    tmpPtr = vv;
    for (i = 0; i < bpp; i++) {
        int j;
        for(j = 0; j < 3; j++) {
            *(tmpPtr++) = 0;
        }
    }
    
    for (it = 0; it < iter; it++) { // iteration
        double ax, ay, az;
        int nit;
        for (nit = 0; nit < nprocs; nit++) { // do nprocs times of comunications between two adjacent processors
            int j;
            for (j = 0; j < bpp; j++) { // for every body in horse
                int k;
                for (k = 0; k < bpp; k++) { // for every body in local
                    double sx = s[k * 3] - ss[j * 3];
                    double sy = s[k * 3 + 1] - ss[j * 3 + 1];
                    double sz = s[k * 3 + 2] - ss[j * 3 + 2];
                    double sd = sqrt(sx * sx + sy * sy + sz * sz);
                    if (sd == 0.0){
                        continue;
                    }
                    sd = sd * sd * sd;
                    ax = G * m[k] * sx / sd;
                    ay = G * m[k] * sy / sd;
                    az = G * m[k] * sz / sd;
                    vv[j * 3] += ax * timestep;
                    vv[j * 3 + 1] += ay * timestep;
                    vv[j * 3 + 2] += az * timestep;
                }
            }
            int next = myrank + 1;
            int last = myrank - 1;
            
            if (next == nprocs) { // last processor
                next = 0;
            }
            if (last < 0) { // first processor
                last = nprocs - 1;
            }
            if (myrank % 2 == 0) {
                MPI_Send(ss, bpp * 3, MPI_DOUBLE, next, TAG3, MPI_COMM_WORLD); // !! tags
                MPI_Send(vv, bpp * 3, MPI_DOUBLE, next, TAG4, MPI_COMM_WORLD); // !! the first parameter
                MPI_Recv(ss, bpp * 3, MPI_DOUBLE, last, TAG5, MPI_COMM_WORLD, &status);
                MPI_Recv(vv, bpp * 3, MPI_DOUBLE, last, TAG6, MPI_COMM_WORLD, &status);
            }
            else {
                double *tmps, *tmpv;
                tmps = (double*)malloc(sizeof(double) * bpp * 3); // allocate a new block of memory to store incmoing date
                tmpv = (double*)malloc(sizeof(double) * bpp * 3);
                
                MPI_Recv(tmps, bpp * 3, MPI_DOUBLE, last, TAG3, MPI_COMM_WORLD,&status); // only passing double pointers
                MPI_Recv(tmpv, bpp * 3, MPI_DOUBLE, last, TAG4, MPI_COMM_WORLD,&status);
                MPI_Send(ss, bpp * 3, MPI_DOUBLE, next, TAG5, MPI_COMM_WORLD);
                MPI_Send(vv, bpp * 3, MPI_DOUBLE, next, TAG6, MPI_COMM_WORLD);
                double *adds = ss;
                double *addv = vv;
                ss = tmps; // move pointers
                vv = tmpv;
                free(adds); // free unneeded memory
                free(addv);
            }
        }
        int j;
        for (j = 0; j < bpp; j++) { // update position coordinates
            ss[j * 3] += vv[j * 3] * timestep;
            ss[j * 3 + 1] += vv[j * 3 + 1] * timestep;
            ss[j * 3 + 2] += vv[j * 3 + 2] * timestep;
            s[j * 3] = ss[j * 3];
            s[j * 3 + 1] = ss[j * 3 + 1];
            s[j * 3 + 2] = ss[j * 3 + 2];
        }
    }
    
    
    // This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
    // in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
    // sure you dont add extra debugging statements in stderr.
    double *newS, *newV, *newM;
    if (myrank == 0) {
        newS = (double *)malloc(sizeof(double) * n * 3);// allocate memory to receive incomming data
        newV = (double *)malloc(sizeof(double) * n * 3);
        newM = (double *)malloc(sizeof(double) * n);
    }
    MPI_Gather(ss, bpp * 3, MPI_DOUBLE, newS + (bpp * 3 * myrank), bpp * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(vv, bpp * 3, MPI_DOUBLE, newV + (bpp * 3 * myrank), bpp * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(m, bpp, MPI_DOUBLE, newM + (bpp * myrank), bpp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (myrank == 0) {
        for (i = 0; i < n; i++) {
            fprintf(stderr, OUTPUT_BODY, newS[i*3], newS[i*3+1], newS[i*3+2], newV[i*3], newV[i*3+1], newV[i*3+2], newM[i]);
        }
        free(newS);
        free(newV);
        free(newM);
    }
    free(ss);
    free(vv);
}
