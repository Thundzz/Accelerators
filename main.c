#include "particule.h"

#include <stdio.h>
#include <stdlib.h>
#define NBITER 1


/* MDTSC */

int main(int argc, char ** argv)
{
	if(argc != 2){
		fprintf(stderr, "Enter the number of particles\n");
		exit(EXIT_FAILURE);
	}
	int NBPAR = atoi(argv[1]);
	double dt = 200.0;
	pset *s = pset_alloc(NBPAR);
	pset_init_orbit(s);

	/*pset_print(s);*/

	FILE * fichier =fopen("datafile", "w+");
	/*fprintf(fichier, "#particule X Y Z\n");*/
	for (int i = 0; i < NBITER ; ++i)
	{
		pset_step(s, dt);
		/*for (int j = 0; j < NBPAR; ++j)
		{
			fprintf(fichier, 
			"%d %g %g %g\n",
			j, s->pos[j], s->pos[j+NBPAR], s->pos[j+2*NBPAR]);
		}
		if(i!= NBITER -1)
			fprintf(fichier, "\n\n");*/
	}

	/*pset_print(s);*/
	

	fclose(fichier);
	pset_free(s);
	return 0;
}