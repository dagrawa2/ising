//      ising.c
// Per compilar:
// gcc -o ising ising.c  -lm 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>

#include "memory.h"
#include "fileio.h"
#include "statistics.h"
#include "random_num.h"


/** Function that shows by screen the help menu. */
void help()
{
    fprintf(stderr,"===============================================\n");
    fprintf(stderr,"Ising model simulation on a hypercubic lattice:\n");
    fprintf(stderr," Usage:\n");
    fprintf(stderr,"  ising [ -d 2 ] [ -L 20 ] [ -J 1 ] [ -B 0.0 ] [ -T 2.0 ] [ [-m | --nmeas ] 1 ] [ [ -c | --nmcs ] 10000000 ] \n");
    fprintf(stderr," Help:\n");
    fprintf(stderr,"  -d [2]: Lattice dimension\n");
    fprintf(stderr,"  -L [20]: Lattice side length\n");
    fprintf(stderr,"  -J [1]: 1 for ferromagnetic, -1 for antiferromagnetic\n");
    fprintf(stderr,"  -B [0.0]: External magnetic field strength\n");
    fprintf(stderr,"  -T [2.0]: Temperature of the simulation\n");
    fprintf(stderr,"  --nmeas [1]: Number of MonteCarlo samples between two measures \n");
    fprintf(stderr,"  --nmcs [10000000]: Total number of MonteCarlo samples\n");
    fprintf(stderr,"  -s [107]: Random Number generator seed\n");
    fprintf(stderr,"  --ieq [1000]: Extra iterations to reach equilibrium\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"===============================================\n");
    fprintf(stderr," Example:\n");
    fprintf(stderr,"  ising -d 4 -L 20 -J 1 -B 0.1 -T 4.0 --nmeas 20 --nmcs 10000000 --ieq 5000\n");
    return;
}

/** Function that analyzes the arguments of the program and saves them 
 * d: dimension of the lattice
 * L: nº of spins at one side of the lattice
 * J: ferromagnetic (1) vs antiferromagnetic (-1)
 * B: External magnetic field strength
 * T: Temperature of the simulation
 * nmeas: Number of MCS between two measures
 * nmcs: Total number of MCS in our simulation
*/
int parseinput(int argc, char **argv,int *d, int *L, int *J, double *B, double *T, \
               int *nmeas,unsigned long int *nmcs, \
               unsigned long int *seed,int *ieq)
{
    int c;
    int itmp;
    unsigned long int ulitmp;
    double dtmp;
    while (1)
    {
        static struct option long_options[] =
        {
            {"d",  required_argument, 0, 'D'},
            {"L",  required_argument, 0, 'L'},
            {"J",  required_argument, 0, 'J'},
            {"B",    required_argument, 0, 'B'},
            {"T",    required_argument, 0, 'T'},
            {"nmeas",    required_argument, 0, 'm'},
            {"nmcs",    required_argument, 0, 'c'},
            {"S", required_argument, 0, 's'},
            {"ieq", required_argument, 0, 'e'},
            {"help", 0, 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long (argc, argv, "hpD:d:l:L:J:j:B:b:T:t:m:c:s:y:", long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c)
        {
            case 'D':
            case 'd':
                if(sscanf(optarg, "%d", &itmp) != 1)
                    fprintf(stderr,"Error parsing -d option\n");
                else
                    *d=itmp;
                break;
            case 'L':
            case 'l':
                if(sscanf(optarg, "%d", &itmp) != 1)
                    fprintf(stderr,"Error parsing -L option\n");
                else
                    *L=itmp;
                break;
            case 'J':
            case 'j':
                if(sscanf(optarg, "%d", &itmp) != 1)
                    fprintf(stderr,"Error parsing -J option\n");
                else if(itmp!=-1 && itmp!=1)
                    fprintf(stderr,"-J option should be -1 or 1\n");
                else
                    *J=itmp;
                break;
            case 'B':
            case 'b':
                if(sscanf(optarg, "%lf", &dtmp) != 1)
                    fprintf(stderr,"Error parsing -B option\n");
                else
                    *B=dtmp;
                break;
            case 'T':
            case 't':
                if(sscanf(optarg, "%lf", &dtmp) != 1)
                    fprintf(stderr,"Error parsing -T option\n");
                else
                    *T=dtmp;
                break;
            case 'm':
                if(sscanf(optarg, "%i", &itmp) != 1 )
                    fprintf(stderr,"Error parsing -m option\n");
                else
                    *nmeas=itmp;
                break;
            case 'c':
                if( sscanf(optarg, "%lud", &ulitmp) != 1 )
                    fprintf(stderr,"Error parsing --nmcs option\n");
                else
                    *nmcs=ulitmp;
                break;
            case 's':
                if( sscanf(optarg, "%lui", &ulitmp) != 1 )
                    fprintf(stderr,"Error parsing -s option\n");
                else
                    *seed=ulitmp;
                break;
            case 'e':
                if(sscanf(optarg, "%i", &itmp) != 1 )
                    fprintf(stderr,"Error parsing --ieq option\n");
                else
                    *ieq=itmp;
                break;
            case 'h':
                help();
                exit(-1);
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                abort();
        }
    }
    return 0;
}


/** Function that gives us the index of a spin given its coordinates */
int xyz2i(const int d, const int L, int *coord)
{
    int i;
    int index=0;
    for (i=0;i<d;++i)
    {
        index+=pow(L,d-i-1)*coord[i];
    }
    return index;
}


/** Function that gives the coordinates for each spin given its number, d and L */
void i2xyz(const int d, const int L, int index, int *coord)
{
    int i;
    for (i=0;i<d;++i)
    {
        coord[i]=index/(int)pow(L,d-i-1);
        index=index-(int)pow(L,d-i-1)*coord[i];
        index=index%((int)pow(L,d-i-1));
    }
    return;
}


/** Function that initializes the lattice */
void InitHypercubicLattice(const int d,const int L,ListIntVector veins)
{
    const int N=pow(L,d); //Number of spins
    int i,j;
    int coord[d];
    IntVector vei;
    for (i=0;i<N;++i)
    {
        veins[i]=CreateIntVector(2*d+1); // In a hypercubic lattice, each spin has 2*d neighbors
        vei=veins[i];
        i2xyz(d,L,i,coord); // Compute the coordinates of this spin
        for (j=0;j<d;++j)
        {
            coord[j]=(coord[j]+1)%L; // Sum 1 to the coordinate j.
            vei[2*j]=xyz2i(d,L,coord); // Neighbor to the index corresponding to this coordinate. 
            coord[j]=(coord[j]-2+L)%L; // Minus 2
            vei[2*j+1]=xyz2i(d,L,coord); // Neighbor
            coord[j]=(coord[j]+1)%L; // Sum 1 (return to the initial position).
        }
        vei[2*d] = N;
    }
    veins[N]=CreateIntVector(N);
    vei=veins[N];
    for (i=0; i<N; ++i)
    {
        vei[i] = i;
    }
    return;
}


/** Create the bonds between neighbors */
void InitBonds(int N, ListIntVector veins, ListIntVector bonds)
{
    int i;
    IntVector vei;
    
    for (i=0;i<=N;++i)
    {
        vei=veins[i];
        bonds[i]=CreateIntVector(vei[-1]);
        ZeroIntVector(bonds[i]);
    }
}


/** Function that initializes the spins in vector v to -1 or +1. */
void RandomVector(IntVector v)
{
    int i;
    for (i=0;i<v[-1];++i)
    {
        v[i]=2*urand(0,2)-1;
    }
    return;
}


/** Function that flips an entire cluster for Wolff algorithm */
void FlipCluster(IntVector spins, IntVector cluster)
{
    int N = cluster[-1]-1;
    int i;
    int j = 1-cluster[N];
    for (i=0;i<spins[-1];++i)
    {
        if (cluster[i] == j) //Pertany al cluster a girar
            spins[i]*=-1;
    }
}


/** Function that gives us if the bond is occupied (1) or not occupied (0)*/
int OccupyBondProb(int Si, int Sj, double expo, int J)
{
    int delta=abs((Si+J*Sj)/2.0);
    if (delta==0 || drand(0.0,1.0) > (1.0-expo)) 
        return 0;
    else
    {
        return 1;
    }
}


/** Function that sets a bond between two spins at the value (usually 1)*/
void SearchBond(IntVector vei, IntVector bond, int spin, int value)
{
    int j;
    for (j=0;j<vei[-1];++j)
    {
        if (spin == vei[j])
            bond[j]=value;
    }
}


/** Function that, for each bond around the spin "around", compute if it is occupied or not*/
void BondsAround(IntVector spins, IntVector cluster, ListIntVector veins, ListIntVector bonds, int around, double expo_B, double expo_J, int J)
{
    int i;
    
    IntVector vei=veins[around]; 
    IntVector bond=bonds[around];
    int N = spins[-1]-1;

    double expo;

    for (i=0;i<vei[-1];++i)
    {
        if (around==N || i==N)
            expo = expo_B;
        else
            expo = expo_J;
        if (bond[i]==0 && (cluster[vei[i]] == 1 || OccupyBondProb(spins[around],spins[vei[i]],expo,J) == 1)) 
        //Ocupem els bonds amb els veins?
        {
            bond[i] = 1;
            //Busquem ara el bond invers per tal de ficar-lo també a 1
            SearchBond(veins[vei[i]],bonds[vei[i]],around,1);
        }
    }
}


/** This recursive function calls itself each time it finds a new spin inside the cluster
 * First, we put the spin inside the cluster
 * Then, we search for occupied bonds between it and its neighbors
 * Finally, for each occupied bond, we move to the corresponding neighbor, EXCEPT when we already visited it.  
 * */

void ClusteringNeigh(IntVector spins, IntVector cluster, ListIntVector veins, ListIntVector bonds, double expo_B, double expo_J, int spin, int J)
{
    int i;
    IntVector vei,bond;
    
    vei=veins[spin];
    bond=bonds[spin];
    
    cluster[spin]=1; //The spin is added to the cluster
    BondsAround(spins,cluster,veins,bonds,spin,expo_B,expo_J,J); //Search for occupied bonds of this spin
    
    for (i=0;i<vei[-1];++i)
    {
        if (bond[i] == 1 && cluster[vei[i]] == 0) //For each bond occupied and unvisited neighbors...
        {
            ClusteringNeigh(spins,cluster,veins,bonds,expo_B,expo_J,vei[i],J);
        }
    }
}


/** Function that perform "steps" clusterings. */
void WolffMetropolis(IntVector spins, ListIntVector veins, ListIntVector bonds, double expo_B, double expo_J, int steps, int J)
{
    int initial;
    int spinscluster;
    int stp;
    int N=spins[-1]-1;
    IntVector cluster=CreateIntVector(N+1);
    int totalspins=0,i;

    for (stp=0;stp<steps;++stp)
    {
        for (totalspins=0;totalspins<N;) //We flip, at least, N spins
        {
            ZeroListIntVector(bonds,N+1); //Sets all bonds to 0 (unoccupied)
            ZeroIntVector(cluster); //Sets the cluster values to 0 (no spins in the cluster)
            initial=urand(0,spins[-1]); //Random spin from which we begin the clustering
            ClusteringNeigh(spins,cluster,veins,bonds,expo_B,expo_J,initial,J); //Create the cluster
            spinscluster=-spins[initial]; //Value of the future spins of the cluster
            for (i=0;i<N;++i) totalspins+=cluster[i]; //Compute the spins in the cluster
            FlipCluster(spins,cluster); //Flip the cluster
        }
    }
    DestroyIntVector(cluster);
}


int main(int argc, char **argv)
{
    int d=2,L=20,nmeas=1; //Dimension of the lattice, size of one side, and Nmeas by default
    unsigned long int nmcs=10000000; //Number of MonteCarlo Samples by default
    int J=1; //+-1 for ferromagnetic vs antiferromagnetic
    double B=0.0; //Default external magnetic field
    double T=2.0; //Temperature by default
    double expo_J; //1-probability
    double expo_B; //1-probability for ghost spin
    int N; //Number of spins
    unsigned long int seed=107; //Seed by default
    int ieq=1000; //Number of MCS, by default, that we ignore by belonging to the termalization
    int i, index;
    int nprog;  // number of MC iterations between printing percentage complete

    char dirname[50];
    char filename1[50];
    char filename2[50];
    
    parseinput(argc,argv,&d,&L,&J,&B,&T,&nmeas,&nmcs,&seed,&ieq); //Function that read the inputs
    printf("Running with parameters: d: %d\t L: %d\t J: %d\t B: %f\t T: %f\t nmeas: %d \t ieq: %d \t nmcs: %ld\n",d,L,J,B,T,nmeas,ieq,nmcs); //Print the inputs

    snprintf(dirname, 50, "output_B%1.3f_T%1.2f", B, T);
    strcpy(filename1, dirname);
    strcpy(filename2, dirname);

    makedir(dirname);
    FILE *json_fitx=Fopen(strcat(filename1, "/args.json"), "w"); // record cmd args in JSON
    fprintf(json_fitx, "{\n  \"d\": %d,\n  \"L\": %d,\n  \"J\": %d,\n  \"B\": %f,\n  \"T\": %f,\n  \"nmeas\": %d,\n  \"nmcs\": %ld,\n  \"s\": %ld,\n  \"ieq\": %d\n}", d, L, J, B, T, nmeas, nmcs, seed, ieq);
    Fclose(json_fitx);

    FILE *states_fitx=Fopen(strcat(filename2, "/states.txt"), "w"); // file to record states

    N=pow(L,d); //Compute N for a square lattice
    expo_J=exp(-2.0/T);
    expo_B=exp(-2.0*B/T);
    InitRNG(seed);
    IntVector spins=CreateIntVector(N+1); //Vector of spins
    ListIntVector veins = CreateListIntVector(N+1); //Vector of neighbors
    ListIntVector bonds; //Bonds between spins for Wolff algorithm
    InitHypercubicLattice(d,L,veins); // Initialize the geometry.
    bonds=CreateListIntVector(N+1);
    InitBonds(N,veins,bonds); //Create the bonds
    RandomVector(spins); // Random initialization of spins
    spins[N] = 1;
    nprog = (nmcs/nmeas)/10;
    for (i=0;i<nmcs/nmeas;++i)
    {
        WolffMetropolis(spins,veins,bonds,expo_B,expo_J,nmeas,J); //Metropolis for "nmcs" MCS taking measures each "nmeas"
        for (index=0; index<N; ++index)
        {
            fprintf(states_fitx, "%d", (spins[index]+1)/2); //write state to file
        }
        if ((i+1)%nprog == 0)
            printf(". . . %d%%\n", 10*(i+1)/nprog); //print progress
    }
    Fclose(states_fitx);
    printf("Done!\n");
    DestroyListIntVector(bonds,N);
    DestroyIntVector(spins);
    DestroyListIntVector(veins,N);
    FreeRNG();
    return 0;
}

