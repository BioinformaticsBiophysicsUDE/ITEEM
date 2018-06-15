// Compile:   g++ ITEEM.cc -std=c++11 -o output.out
// Run:       ./output.out $lambda $delta $tstart $tfinal $step_write            (0<delta<1 ; lambda= integer and equal to 0 for infinite lifespan ; tstart must be 0 for the first run)

#include <iostream>
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
#include <fstream>
#include <string>       // strings
#include <sstream>      // stringstreams
#include <random>
#include <chrono>
#include <math.h>
#include <cmath>

using namespace std;

const int Nmax=100000;
// Nmax: Number of sites
const double fraction=0.5;
// fraction: fraction of sites occupied in the first step
const int Nin=1000;
// Nin: maximum allowed number of strains (Increasing Nin will slow down the program but by increasing the mutation_rate,
// inevitably the number of different species will increase. So should be adjusted accordingly manually.
const double mutation_rate=0.001;
// Mutation rate
const double mutation_width=0.02;
// Sigma of the gaussian distribution for mutation
long int lambda; 
// Average of the attributed lifespan (age of natural death); if set to 0 then there is no attributed lifespan (infinite lifespan)
double delta;
// Trade-off parameter (greater than 0 and less than 1)
double m_gamma;
// Trade-off exponent


float competition[Nin+1][Nin+1];
// Matrix for the competition; 0.5*(1+competition[i][j]) is the interaction matrix introduced in the manuscript = the probability that i can win in the competition with j

int tstart,tfinal,step_write;
// tstart=the generation that the simulation starts from
// tfinal=the generation that the simulation stops at
// Every step_write generations the results are written in a file

int Type[Nmax+1],Age[Nmax+1],lifespan[Nmax+1];
// Type shows that individual number i belongs to which spieaces
// Age stores the age of individuals
// lifespan stores the attributed lifespan of individuals

int N;
// Number of strains at each time step

int abundancy[Nin+1];
//abundancy of strains; 
int first_parent[Nin+1],generation[Nin+1];
//first_parent stores the first parent of strains, it is relevant when we strat from a realization with several strains. In our simulations (speciation from 1 common ancestor) it is 1 for all strains.
//generation stores the position of strains in the geneaological tree
int mutation_time[Nin+1];
//mutation_time stores the generation at which the strains emerge, 
float reproduction_rate[Nin+1],sum_ability[Nin+1];
//reproduction rate and sum of ability of strains;  (competitive ability in the manuscript is the average of sum_ability)
unsigned int generation_history[Nin+1][Nin+1],real_ID[Nin+1],parent[Nin+1],newborn;
//generation_history stores the list of all ancestors of strains
//real_ID stores the order of emergence of strains
//parent stores parent of strains using the real_ID
//newborn stores the real_ID of the last emerged strain
unsigned int list_id[Nin+1];
// To speed-up the simulation list_id stores the id of arrays that contain strains (not extincted)
double reproduce_counter=0,invade_counter=0,death_counter=0;
// They store the number of reproductions, invasions and deaths in each generation

int tt;
//time, each time_step is consisted of Nmax operation of choosing two individual
ofstream myWriteFile_t;

///////////////////////////////////////////////////

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
uniform_int_distribution<int> distribution(1,Nmax);

///////////////////////////////////////////////////

void MutateCompetition(int spnode, int nspnode);
void initialize_configuration();
void initialize_competition();
void realize_initial_configuration_one_type();

///////////////////////////////////////////////////  

double gaussian(double mean, double sigma) {
  // generates gaussian distribution according to Box-Muller algorithm.
  double x1,x2;
  x1=((double)(rand())+1.)/((double)(RAND_MAX)+1.);
  x2=((double)(rand())+1.)/((double)(RAND_MAX)+1.);
  return(mean+sigma*sqrt(-2.*log(x1))*cos(2.*M_PI*x2));
}

///////////////////////////////////////////////////
// Updates reproduction rate after changes in the competition matrix
void update_reproduction() {
  double comp_ab;
  if (N>1) {
    for(int kd=1; kd<=N; kd++) {
      int k=list_id[kd];
      comp_ab=0.5+sum_ability[k]/(2.*(N-1)); //competitive ability
      reproduction_rate[k]=pow((1.-pow(comp_ab,m_gamma)),(1./m_gamma)); // Eq. 3 in the appendix
    }
  } else {
    reproduction_rate[list_id[1]]=1.;
  }
}

///////////////////////////////////////////////////  
// Removes the competition rows and columns and all of the other attributes of a given strain. It is called when the abundancy of one strain goes to 0
void remove(int ss) {
  for(int kd=1; kd<=N; kd++) {
    int k=list_id[kd];
    sum_ability[k]=sum_ability[k]-competition[k][ss];
    competition[k][ss]=0;
    competition[ss][k]=0;
  }
  sum_ability[ss]=0;
  reproduction_rate[ss]=0;
  real_ID[ss]=0;
  first_parent[ss]=0;
  parent[ss]=0;
  for (int k=1; k<=generation[ss]; k++) {
    generation_history[ss][k] = 0;
  }
  generation[ss]=0;
  if (tt>mutation_time[ss]+10) {
    myWriteFile_t << mutation_time[ss] << "\t" << tt << "\n";
  }
  mutation_time[ss]=0;
  update_reproduction();
}

///////////////////////////////////////////////////
// Deletes an individual (with id ii and type (strain) sp) and makes sure if it is extinct, the corresponding strain is eliminated
void eliminate(int sp, int ii) {
  Type[ii]=0;
  Age[ii]=0;
  lifespan[ii]=0;
  abundancy[sp]--;
  if(abundancy[sp]==0) {
    N--;
    if (N>0) {
      int id=1;
      for (int i=1; i<=Nin; i++) {
	if (abundancy[i]!=0) {
	  list_id[id]=i;
	  id++;
	}
      }
      remove(sp);
    } else {
      initialize_configuration();
      initialize_competition();
      realize_initial_configuration_one_type();
    }
  }
}

///////////////////////////////////////////////////
// Finds an ID ready to be used as new species
int subNew() {
  int nsp=1;
  for(nsp=1; nsp<=Nin; nsp++) {
    if(abundancy[nsp]==0) {
      break;
    }
  }
  if(nsp>Nin) {
    cout << "strange happened!\n";
    exit(0);
    // If it happens it means that the maximum number of species (Nin) is not sufficient for the chosen mutation rate
  }
  return(nsp);
}

///////////////////////////////////////////////////
// For a selected id (ii) which is going to be occupied with a mutant of a previously existed type (sp) finds and generates its attributes
void mutate(int sp, int ii) {
  std::poisson_distribution<int> dist_poisson(lambda);

  int nsp=subNew();
  abundancy[nsp]++;
  Type[ii]=nsp;
  Age[ii]=1;
  lifespan[ii]=dist_poisson(generator);
  parent[nsp]=real_ID[sp];
  first_parent[nsp]=first_parent[sp];
  newborn++;
  real_ID[nsp]=newborn;
  generation[nsp]=generation[sp]+1;
  for (int k=0; k<=generation[sp]; k++) {
    generation_history[nsp][k]=generation_history[sp][k];
  }
  generation_history[nsp][generation[nsp]]=real_ID[nsp];
  mutation_time[nsp]=tt;
  N++;
  list_id[N]=nsp;
  MutateCompetition(sp,nsp);
}

///////////////////////////////////////////////////
// Puts one individual of tyepe sp in the cite ii (simple reproduction)
void reproduce(int sp, int ii) {
  std::poisson_distribution<int> dist_poisson(lambda);

    abundancy[sp]++;
    Type[ii]=sp;
    Age[ii]=1;
    lifespan[ii]=dist_poisson(generator);
}

///////////////////////////////////////////////////
// Invasion step
// Selects one individual randomely then it reproduces based on its reproduction rate
// Selects another individual randomely then calculates the probability to win if it has been occupied by another individual
// calculates the mutation probability for each invasion

void invasion() {
  int ii,jj;
  int spec1,spec2; 
  double prob;
  ii = distribution(generator);
  spec1=Type[ii];
  if (((double)rand()/(double)RAND_MAX)<reproduction_rate[spec1]) {
    jj=distribution(generator);
    while(ii==jj) 
      jj=distribution(generator);
    spec2=Type[jj];
    if(spec2==0) {
      reproduce_counter++;
      if (((double)rand()/(double)RAND_MAX)<mutation_rate) {
	mutate(spec1,jj);
      } else {
	reproduce(spec1,jj);
	
      }    
    } else {
      prob=0.5*(1+competition[spec1][spec2]);
      if (((double)rand()/(double)RAND_MAX)<prob) { 
	eliminate(spec2,jj);
	invade_counter++;
	if (((double)rand()/(double)RAND_MAX)<mutation_rate) {
	  mutate(spec1,jj);
	} else {
	  reproduce(spec1,jj);
	}
      }
      
    }
  }
}

//////////////////////////////////////////////////////
// Output in file
void WriteAge(int iitime) {
  stringstream file_name;
  file_name << "./" << lambda << "/" << delta << "/" << iitime << "_Age_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ofstream myWriteFile;
  myWriteFile.open(file_name.str().c_str());
  for (int i=1; i<=Nmax; i++) {
    myWriteFile << Type[i] << "\t" << Age[i] << "\t" << lifespan[i] << "\n";
  }
  myWriteFile.close();
}

///////////////////////////////////////////////////
// Output in file
void WriteCompetition(int iitime) {
  stringstream file_name;
  file_name << "./" << lambda << "/" << delta << "/" << iitime << "_competition_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ofstream myWriteFile;
  myWriteFile.open(file_name.str().c_str());
  myWriteFile << "i\treal_ID\tparent\t1stPar\tgenera\tmutation\tabundancy\tbirth\tdeath\tability\n";
  myWriteFile << "\t\t\t\t\t\t\t\t\t\t\t";
  for (int id=1; id<=N; id++) {
    myWriteFile << list_id[id]  << "\t";
  }
  myWriteFile << "\n";
  for (int id=1; id<=N; id++) {
    int i=list_id[id];
    myWriteFile << i << "\t" << real_ID[i] << "\t" << parent[i] << "\t" << first_parent[i] << "\t" << generation[i] << "\t" << mutation_time[i] << "\t" << abundancy[i] << "\t" << reproduction_rate[i] << "\t" << sum_ability[i] << "\t";
    for (int jd=1; jd<=N; jd++) {
      myWriteFile << competition[i][list_id[jd]] << "\t";
    }
    myWriteFile << "\n";
  }
  myWriteFile.close();
} 

///////////////////////////////////////////////////
// Output in file
void WriteGeneration(int iitime) {
  stringstream file_name;
  file_name << "./" << lambda << "/" << delta << "/" << iitime << "_generation_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ofstream myWriteFile;
  myWriteFile.open(file_name.str().c_str());
  for (int id=1; id<=N; id++) {
    int i=list_id[id];
    myWriteFile << i << "\t";
    for (int k=1; k<=generation[i]; k++) {
      myWriteFile << generation_history[i][k] << "\t";
    }
    myWriteFile << "\n";
  }
  myWriteFile.close();
} 

///////////////////////////////////////////////////
// Output in file
void WriteList(int iitime) {
  stringstream file_name;
  file_name << "./" << lambda << "/" << delta << "/" << iitime << "_List_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ofstream myWriteFile;
  myWriteFile.open(file_name.str().c_str());
  for (int id=1; id<=N; id++) {
    myWriteFile << id << "\t" << list_id[id] << "\n";
  }
  myWriteFile.close();
} 

///////////////////////////////////////////////////
// Reads initial configuration to restart a simulation from certain time
void ReadConfiguration(int iitime) {

  stringstream file_name2;
  file_name2 << "./" << lambda << "/" << delta << "/" << tstart << "_competition_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ifstream myReadFile2;
  myReadFile2.open(file_name2.str().c_str());
  
  string line;
  int dummy;
  int jj,ii;
  int species[Nin+1];
  jj=1;
  getline(myReadFile2,line);
  getline(myReadFile2,line);
  istringstream iss( line );
  while( iss >> ii ) {
    species[jj]=ii;
    jj++;
  }
  N=jj-1;
  double d;
  while(getline(myReadFile2,line)) {
    istringstream iss(line);
    iss >> ii;
    iss >> real_ID[ii];
    iss >> parent[ii];
    iss >> first_parent[ii];
    iss >> generation[ii];
    iss >> mutation_time[ii];
    iss >> abundancy[ii];
    iss >> reproduction_rate[ii];
    iss >> sum_ability[ii];
    jj=1;
    while( iss >> d ) {
      competition[ii][species[jj]]=d;
      jj++;
    }
  }

  stringstream file_name3;
  file_name3 << "./" << lambda << "/" << delta << "/" << tstart << "_Age_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ifstream myReadFile3;
  myReadFile3.open(file_name3.str().c_str());


  jj=1;
  while(getline(myReadFile3,line)) {
    istringstream iss(line);
    iss >> Type[jj];
    iss >> Age[jj];
    iss >> lifespan[jj];
    jj++;
  }
  
  stringstream file_name4;
  file_name4 << "./" << lambda << "/" << delta << "/" << tstart << "_generation_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ifstream myReadFile4;
  myReadFile4.open(file_name4.str().c_str());
  while(getline(myReadFile4,line)) {
    istringstream iss(line);
    iss >> ii;
    for (int k=1; k<=generation[ii]; k++) {
      iss >> generation_history[ii][k];
    }
  }
  int id;
  stringstream file_name5;
  file_name5 << "./" << lambda << "/" << delta << "/" << tstart << "_List_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ifstream myReadFile5;
  myReadFile5.open(file_name5.str().c_str());
  while(getline(myReadFile5,line)) {
    istringstream iss(line);
    iss >> id;
    iss >> list_id[id];
  }
}

///////////////////////////////////////////////////
// initializes the variables; sets everything to zero
void initialize_configuration() {
  for (int i=0; i<=Nmax; i++) {
    Type[i]=0;
    Age[i]=0;
    lifespan[i]=0;
  }
  for (int i=0; i<=Nin; i++) {
    abundancy[i]=0;
    real_ID[i]=0;
    list_id[i]=0;
    first_parent[i]=0;
    parent[i]=0;
    generation[i]=0;
    mutation_time[i]=0;
    reproduction_rate[i]=0;
    sum_ability[i]=0;
    for (int j=0; j<=Nin; j++) {
      generation_history[i][j]=0;
    }
  }
  N=0;
  newborn=0;
}  

///////////////////////////////////////////////////
// initializes the competition matrix; sets everything to zero
void initialize_competition() {
  for (int i=1; i<=Nin; i++) {
    for (int j=1; j<=Nin; j++) {
      competition[i][j]=0.0;
    }
  }
}  

///////////////////////////////////////////////////
// Fills the fraction of sites with 1 strain (ancestor of all)
void realize_initial_configuration_one_type() {
  std::poisson_distribution<int> dist_poisson(lambda);
  
  N=1;
  list_id[1]=1;
  parent[1]=0;
  first_parent[1]=1;
  for (int i=1; i<=int(Nmax*fraction); i++) {
    Type[i]=1;
    Age[i]=1;
    lifespan[i]=dist_poisson(generator);
  }
  generation[1]=1;
  generation_history[1][1]=1;
  abundancy[1]=Nmax/2;
  MutateCompetition(1,1);
  newborn++;
  real_ID[1]=newborn;
  sum_ability[1]=0;
}  

///////////////////////////////////////////////////
// produces the interaction trait (row and column) of the new emerged mutant
void MutateCompetition(int spnode, int nspnode) {
  for(int jd=1;jd<=N;jd++) {
    int j=list_id[jd];
    if (j!=nspnode) {
      competition[nspnode][j]=competition[spnode][j]+gaussian(0,mutation_width);
      if (competition[nspnode][j]>1) {
	competition[nspnode][j]=2-competition[nspnode][j];
      } else if (competition[nspnode][j]<-1) {
	competition[nspnode][j]=-2-competition[nspnode][j];
      }
      competition[j][nspnode]=-competition[nspnode][j];
      sum_ability[j]=sum_ability[j]+competition[j][nspnode];
      sum_ability[nspnode]=sum_ability[nspnode]+competition[nspnode][j];
    }
  }
  competition[nspnode][nspnode]=0;
  update_reproduction();
}
  

///////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  
  lambda = atof(argv[1]); //reads the lambda from the input
  delta = atof(argv[2]);  //reads the delta from the input
    
  tstart = atof(argv[3]); //reads the value of the first step from the input
  tfinal = atof(argv[4]); //reads the value of the last step from the input
  step_write = atof(argv[5]); //reads the value of the step_write from the input

  m_gamma = (-log(1-delta))/0.693147181;   // -log ( 1 - delta ) / log(2) = log2 ( 1 - delta ) ; calculates trade-off exponent from trade-off parameter
   
  int cell;
  int seed2 = time(NULL);
  srand (seed2);
  
  // file that stores the first and the last generations that a type (strain) exists
  stringstream file_name_t;
  file_name_t << "./" << lambda << "/" << delta << "/" << tfinal << "_Duration_" << Nmax << "_" << Nin << "_" << mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  myWriteFile_t.open(file_name_t.str().c_str());

  initialize_competition();
  
  //file that keeps track of the reproduction, death and succesful invasions
  stringstream file_name_r;
  file_name_r << "./" << lambda << "/" << delta << "/" << tfinal << "_Rate_" << Nmax << "_" << Nin << "_" << 
  mutation_rate << "_" << mutation_width << "_" << lambda << "_" << delta << ".dat";
  ofstream myWriteFile_r;
  myWriteFile_r.open(file_name_r.str().c_str());
  
  initialize_configuration();

  if (tstart == 0) {
    realize_initial_configuration_one_type();
  } else {
    ReadConfiguration(tstart);
  }  
  
  // itime=counter for outer step (goes from tstart to tfinal); iinternal=counter for inner step (goes from 1 to step_write)
  for(int itime=tstart/step_write+1; itime<=tfinal/step_write; itime++) {
    if (N > 0) {
      reproduce_counter=0;
      invade_counter=0;
      death_counter=0;
      for(int iinternal=1; iinternal<=step_write; iinternal++) {
	tt = (itime-1)*step_write+iinternal;
	for(cell=1; cell<=Nmax; cell++) {
	  invasion();
	}
	for(cell=1; cell<=Nmax; cell++) {
	  if (Type[cell]!=0) {
	    Age[cell]++;
	    if (Age[cell]>lifespan[cell] && lambda > 0) {
	      eliminate(Type[cell],cell);
	      death_counter++;
	    }
	  }
	}
      }
      cout << "generation = " << itime*step_write << "\n";
      WriteCompetition(itime*step_write);
      WriteGeneration(itime*step_write);
      WriteList(itime*step_write);
      myWriteFile_r << tt << "\t" << step_write << "\t" << reproduce_counter << "\t" << invade_counter << "\t" << death_counter << "\n"; 
    }
  }
  WriteAge(tfinal);
  myWriteFile_t.close();
  myWriteFile_r.close();
  return 0;
}