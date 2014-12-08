#include "particule.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>

#define NBITER 100
#define NBPAR 3

static int SEEDED =0;
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

pset * pset_alloc(int nb_par){
	pset * set = (pset *)malloc(sizeof(pset));
	if(set == NULL)
	{
		fprintf(stderr, "Can't allocate memory for the set creation.\n");
		exit(EXIT_FAILURE);
	}
	set->nb = nb_par;
	set->m = (double*)malloc(nb_par * sizeof(double));
	set->pos = (double*)malloc(3* nb_par * sizeof(double));
	set->spd = (double*)malloc(3* nb_par * sizeof(double));
	set->acc = (double*)malloc(3* nb_par * sizeof(double));
	set->force = (double*)malloc(3* nb_par * sizeof(double));
	return set;
}

void pset_free(pset * set){
	free(set->pos);
	free(set->spd);
	free(set->acc);
	free(set->m);
	free(set->force);
	free(set);
}

void pset_copy(pset * origin, pset * dest){
	int nb = origin-> nb;
	int sd = sizeof(double);
	dest->nb = origin->nb;
	memcpy(dest->m, origin->m  , nb*sd);
	memcpy(dest->acc, origin->acc, 3* nb*sd);
	memcpy(dest->spd, origin->spd, 3* nb*sd);
	memcpy(dest->pos, origin->pos, 3* nb*sd);
}

void pset_print(pset * set)
{
	int i;
	int size = set->nb;
	for (i = 0; i < size; ++i)
	{
		printf("#Particule numéro : %d, de masse %g\n", i, set->m[i]);
		printf("\tx:%g y:%g z:%g\n", set->pos[i], set->pos[i+ size], set->pos[i+ 2*size] );
		printf("\tvx:%g vy:%g vz:%g\n",set->spd[i], set->spd[i+ size], set->spd[i+ 2*size]);
		printf("\tax:%g ay:%g az:%g\n",set->acc[i], set->acc[i+ size], set->acc[i+ 2*size]);
	}
}

void seed()
{
	if(!SEEDED)
	{
		unsigned long seed = mix(clock(), time(NULL), getpid());
		srand(seed);
		SEEDED++;
	}
}

void pset_init_rand(pset * s)
{
	seed();
	int i;
	int size = s->nb;
	for (i = 0; i < size; i++)
	{
		s->m[i] = 1.0e10;
		s->pos[i] = MIN_RAND + rand()%(MAX_RAND-MIN_RAND);
		s->pos[i+size] = MIN_RAND + rand()%(MAX_RAND-MIN_RAND);
		s->pos[i+2*size] = MIN_RAND + rand()%(MAX_RAND-MIN_RAND);
		s->spd[i] = 0;
		s->spd[i+size] = 0;
		s->spd[i+2*size] = 0;
		s->acc[i] = 0;
		s->acc[i+size] = 0;
		s->acc[i+2*size] = 0;
	}
}

/** Prend en argument les masses de deux particules, la distance entre
 * ces particules et retourne l'intensité de la force gravitationnelle
 * entre ces deux particules
 */


__device__ void distance(double x1, double y1, double z1 , double x2, double y2, double z2, double *res)
{
	*res = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
}

__device__ void intensity(double m, double d, double * res)
{
	*res = (CONST_GRAV * m / (d*d*d));
}

__global__ void nbody(int* n, double* acc, double* spd, double* pos, double* m)
{
	unsigned int idx = threadIdx.x;

	int j;
	double d, inten1, inten2;
	int size = *n;
	double dt = 500.0;

	if (m[idx] != 0)
	{
		for (j = idx+1; j < size; ++j)
		{
			if(m[j] != 0)
				distance(pos[idx], pos[idx+ size], pos[idx+ 2*size], pos[j],
						 pos[j+ size], pos[idx+ 2*size], &d);
			else 
				d = DBL_MAX;

			intensity(m[j], d, &inten1);
			acc[idx]+= inten1 *(pos[j] - pos[idx]); 
			acc[idx+size]+= inten1 *(pos[j+size] - pos[idx+size]);
			acc[idx+2*size]+= inten1 *(pos[j+2*size] - pos[idx+2*size]);

			intensity(m[idx], d, &inten2);
			acc[j]-= inten2 *(pos[j] - pos[idx]);  
			acc[j+size]-= inten2 *(pos[j+size] - pos[idx+size]);
			acc[j+2*size]-= inten2 *(pos[j+2*size] - pos[idx+2*size]);
		}
	}
	pos[idx]+= dt* spd[idx] + dt*dt/2 * acc[idx];
	pos[idx + size]+= dt* spd[idx+ size] + dt*dt/2 * acc[idx+size];
	pos[idx + 2*size]+= dt* spd[idx+ 2*size] + dt*dt/2 * acc[idx+2*size];

	spd[idx]+= dt* acc[idx];
	spd[idx + size]+= dt* acc[idx+ size];
	spd[idx + 2*size]+= dt* acc[idx+ 2*size];
}


int main()
{
	pset *s = pset_alloc(NBPAR);
	pset_init_rand(s);

	pset_print(s);

	int* nb;
	double* acc, *spd, *pos, *m;
	cudaMalloc((void**)&nb, 1*sizeof(int));
	cudaMalloc((void**)&acc, 3*NBPAR*sizeof(double));
	cudaMalloc((void**)&spd, 3*NBPAR*sizeof(double));
	cudaMalloc((void**)&pos, 3*NBPAR*sizeof(double));
	cudaMalloc((void**)&m, NBPAR*sizeof(double));

	cudaMemcpy(nb, &s->nb, 1*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(acc, s->acc, 3*NBPAR*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(spd, s->spd, 3*NBPAR*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(pos, s->pos, 3*NBPAR*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(m, s->m, NBPAR*sizeof(double), cudaMemcpyHostToDevice);

	int dimBlock = 1;


	FILE * fichier =fopen("datafile", "w+");
	fprintf(fichier, "#particule X Y Z\n");
	for (int i = 0; i < NBITER ; ++i)
	{
		nbody<<<1, dimBlock>>>(nb, acc, spd, pos, m);
		cudaMemcpy(s->pos, pos, 3*NBPAR*sizeof(double), cudaMemcpyDeviceToHost);
		for (int j = 0; j < NBPAR; ++j)
		{
			fprintf(fichier, 
			"%d %g %g %g\n",
			j, s->pos[j], s->pos[j+NBPAR], s->pos[j+2*NBPAR]);
		}
		if(i!= NBITER -1)
			fprintf(fichier, "\n\n");
	}

	pset_print(s);
	

	fclose(fichier);
	pset_free(s);
	cudaFree(nb);
	cudaFree(acc);
	cudaFree(spd);
	cudaFree(pos);
	cudaFree(m);
	return 0;
}