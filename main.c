#include "particule.h"

#include <stdio.h>

#define NBITER 100
#define NBPAR 3


/* MDTSC */

int main()
{
	double dt = 500.0;
	pset *s = pset_alloc(NBPAR);
	pset_init_rand(s);

	pset_print(s);

	FILE * fichier =fopen("datafile", "w+");
	fprintf(fichier, "#particule X Y Z\n");
	for (int i = 0; i < NBITER ; ++i)
	{
		pset_step(s, dt);
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
	return 0;
}