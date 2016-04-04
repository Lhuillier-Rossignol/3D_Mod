#include "lib_3d.h"
#include "lib_mat.h"

void multiplicationVecteur3d(t_point3d *v1, double m[4][4], t_point3d *v2) // v1 = m*v2
{
	int cpt;
	int i, j;
	for (i=0;i<4;i++) {
		cpt=0;
		for (j=0;j<4;j++) {
			cpt =cpt+ m[i][j]*(v2->xyzt[j]);
		}
		v1->xyzt[i]=cpt;
		
	}
}

void multiplicationMatrice3d(double m1[4][4], double m2[4][4], double m3[4][4]) // m1 = m2*m3
{
	int i, j, k;
	
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
	  		for (k=0;k<4;k++){ 
	    			// multiplication revien à l'addition de m1[i][k] et de m2[k][j] avec k allant de 0 colonnes
	    			m1[i][j]=m1[i][j]+m2[i][k]*m3[k][j];
	  		}
		}
      	}
}

void copierMatrice3d(double m1[4][4], double m2[4][4]) // m1 = m2
{
	int i, j;

	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			m1[i][j] = m2[i][j];
		}
	}
}
