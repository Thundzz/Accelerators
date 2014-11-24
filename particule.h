#ifndef PARTICULE_H
#define PARTICULE_H


#define CONST_GRAV 6.6742e-11
#define MIN_RAND -1000
#define MAX_RAND 1000

struct p_set {
	int nb;
	double * m;
	double * acc;
	double * spd;
	double * pos;
	double * force;
	int * globId;
};

typedef struct p_set pset;


/** Alloue un ensemble de nb_par particules.
 **/
pset * pset_alloc(int nb_par);


/** Copie le contenu du pset origin dans le pset dest
 **/
void pset_copy(pset * origin, pset * dest);

/** Désallloue un ensemble de particules.
**/
void pset_free(pset * set);



/** Initialise aléatoirement un set
**/
void pset_init_rand(pset * s);



/** Affiche les infos sur les particules d'un ensemble
**/
void pset_print(pset * s);

/** Réalise un pas d'itération sur tous les éléments de l'ensemble s:
* La position de chaque particule est mise à jour correctement, en 
* supposant que la valeur actuelle de l'accéléréation a déjà été
* calculée. i.e que f_grav a déjà été executée avec tous
* les autres ensembles de particules.
**/
void pset_step(pset * s, double dt);


#endif