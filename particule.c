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
	set->m = malloc(nb_par * sizeof(double));
	set->pos = malloc(3* nb_par * sizeof(double));
	set->spd = malloc(3* nb_par * sizeof(double));
	set->acc = malloc(3* nb_par * sizeof(double));
	set->force = malloc(3* nb_par * sizeof(double));
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

/* Calcule la vitesse de satellisation */
double v_orbit(double mass, double distance)
{
	return sqrt(CONST_GRAV*mass/distance);
}

void pset_init_orbit(pset *s)
{
	seed();
	double dmin= 200, distance;
	int size = s->nb;
	s->pos[0 ] = 0;
	s->pos[0 +size] = 0;
	s->pos[0 +2*size] = 0;
	s->spd[0 ] = 0;
	s->spd[0 +size] = 0;
	s->spd[0 +2*size] = 0;
	s->acc[0 ] = 0;
	s->acc[0 + size] = 0;
	s->acc[0 +2*size] = 0;

	s->m[0] = 1e10;

	for (int i = 1; i < size; ++i)
	{
		distance = dmin*i +  rand()% 50;
		s->m[i] = s->m[0] /20000;
		s->pos[i] = s->pos[0] -distance;
		s->pos[i+size] = 0;
		s->pos[i+2*size] = 0;
		s->spd[i] = 0;
		s->spd[i+size]= v_orbit(s->m[0], distance);
		s->spd[i+2*size]= 0;
		s->acc[i] = 0;
		s->acc[i+size] = 0;
		s->acc[i+2*size] = 0;
	}
}

/** Calcule la distance entre deux particules données par leurs
 * coordonnées cartésiennes.
 */
double distance(double x1, double y1, double z1 , double x2, double y2, double z2)
{
	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
}

/** Prend en argument les masses de deux particules, la distance entre
 * ces particules et retourne l'intensité de la force gravitationnelle
 * entre ces deux particules
 */
double intensity(double m, double d)
{
	return (CONST_GRAV * m / (d*d*d));
}

void update_acc(pset* s)
{
	int i, j;
	double d, inten1, inten2;
	int size = s->nb;
	#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < size; ++i)
	{
		s->acc[i] = 0 ; 
		s->acc[i+size] = 0;
		s->acc[i+2*size] = 0;
	}
	#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < size; ++i)
	{
		for (j = i+1; j < size; ++j)
		{
			d = distance(s->pos[i], s->pos[i+ size], s->pos[i+ 2*size],
						 s->pos[j], s->pos[j+ size], s->pos[j+ 2*size]);


			inten1 = intensity(s->m[j], d);
			s->acc[i]+= inten1 *(s->pos[j] - s->pos[i]); 
			s->acc[i+size]+= inten1 *(s->pos[j+size] - s->pos[i+size]);
			s->acc[i+2*size]+= inten1 *(s->pos[j+2*size] - s->pos[i+2*size]);

			inten2 = intensity(s->m[i], d);
			s->acc[j]-= inten2 *(s->pos[j] - s->pos[i]);  
			s->acc[j+size]-= inten2 *(s->pos[j+size] - s->pos[i+size]);
			s->acc[j+2*size]-= inten2 *(s->pos[j+2*size] - s->pos[i+2*size]);
		}
	}
}

void update_spd(pset * s, double dt){
	int size = s->nb;
	for (int i = 0; i < size; ++i)
	{
		s->spd[i]+= dt* s->acc[i];
		s->spd[i + size]+= dt* s->acc[i+ size];
		s->spd[i + 2*size]+= dt* s->acc[i+ 2*size];
	}
}

void update_pos(pset * s, double dt){
	int size = s->nb;
	for (int i = 0; i < size; ++i)
	{
		s->pos[i]+= dt* s->spd[i] + dt*dt/2 * s->acc[i];
		s->pos[i + size]+= dt* s->spd[i+ size] + dt*dt/2 * s->acc[i+size];
		s->pos[i + 2*size]+= dt* s->spd[i+ 2*size] + dt*dt/2 * s->acc[i+2*size];
	}
}


void pset_step(pset * s, double dt)
{
	update_acc(s);
	update_pos(s, dt);
	update_spd(s, dt);
}