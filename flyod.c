#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <mpi.h>

const int INFINITY = 1000000;
void Read_matrix(int loca_lmat[], int n, int myrank, int p,
MPI_Comm comm); //this function reads matrix data from local matrix. void Print_matrix(int loca_lmat[], int n, int myrank, int p,
MPI_Comm comm);
void Floyd(int loca_lmat[], int n, int myrank, int p, MPI_Comm comm); int Owner(int k, int p, int n);
void Copy_row(int loca_lmat[], int n, int p, int ro_wk[], int k);
void Printrow(int loca_lmat[], int n, int myrank, int i); //prints row of the matrix int main(int argc, char* argv[]) {
int n;
int* loca_lmat; MPI_Comm comm; int p, myrank;
MPI_Init(&argc, &argv);
comm = MPI_COMM_WORLD; MPI_Comm_size(comm, &p); MPI_Comm_rank(comm, &myrank); if (myrank == 0) {
printf("How many vertices:-\n"); scanf("%d", &n);
}
MPI_Bcast(&n, 1, MPI_INT, 0, comm); loca_lmat = malloc(n*n/p*sizeof(int));
if (myrank == 0) printf("Enter the loca_lmatrix\n"); Read_matrix(loca_lmat, n, myrank, p, comm);
if (myrank == 0) printf("We got\n"); Print_matrix(loca_lmat, n, myrank, p, comm); if (myrank == 0) printf("\n"); Floyd(loca_lmat, n, myrank, p, comm);
if (myrank == 0) printf("The solution is:\n");
 
Print_matrix(loca_lmat, n, myrank, p, comm); free(loca_lmat);
MPI_Finalize(); return 0;
} /* main */
void Read_matrix(int loca_lmat[], int n, int myrank, int p, MPI_Comm comm) {
int i, j;
int* tempmat = NULL; if (myrank == 0) {
tempmat = malloc(n*n*sizeof(int)); for (i = 0; i < n; i++)
for (j = 0; j < n; j++)
scanf("%d", &tempmat[i*n+j]); MPI_Scatter(tempmat, n*n/p, MPI_INT,
loca_lmat, n*n/p, MPI_INT, 0, comm); free(tempmat);
} else {
MPI_Scatter(tempmat, n*n/p, MPI_INT, loca_lmat, n*n/p, MPI_INT, 0, comm);
}
}
void Printrow(int loca_lmat[], int n, int myrank, int i){ char charint[100];
char charrow[1000]; int j, offset = 0;
for (j = 0; j < n; j++) {
if (loca_lmat[i*n + j] == INFINITY) sprintf(charint, "i ");
else
sprintf(charint, "%d ", loca_lmat[i*n + j]); sprintf(charrow + offset, "%s", charint); offset += strlen(charint);
}
printf("Proc %d > row %d = %s\n", myrank, i, charrow);
} /* Printrow */
 
void Print_matrix(int loca_lmat[], int n, int myrank, int p, //prints the matrix MPI_Comm comm) {
int i, j;
int* tempmat = NULL; if (myrank == 0) {
tempmat = malloc(n*n*sizeof(int)); MPI_Gather(loca_lmat, n*n/p, MPI_INT,
tempmat, n*n/p, MPI_INT, 0, comm); for (i = 0; i < n; i++) {
for (j = 0; j < n; j++)
if (tempmat[i*n+j] == INFINITY) printf("i ");
else
printf("%d ", tempmat[i*n+j]); printf("\n");
}
free(tempmat);
} else {
MPI_Gather(loca_lmat, n*n/p, MPI_INT, tempmat, n*n/p, MPI_INT, 0, comm);
}
}
void Floyd(int loca_lmat[], int n, int myrank, int p, MPI_Comm comm) { //implement Floyd algo int globalk, locali, globalj, temp;
int root;
int* ro_wk = malloc(n*sizeof(int));
for (globalk = 0; globalk < n; globalk++) { root = Owner(globalk, p, n);
if (myrank == root)
Copy_row(loca_lmat, n, p, ro_wk, globalk); MPI_Bcast(ro_wk, n, MPI_INT, root, comm); for (locali = 0; locali < n/p; locali++)
 
for (globalj = 0; globalj < n; globalj++) {
temp = loca_lmat[locali*n + globalk] + ro_wk[globalj]; if (temp < loca_lmat[locali*n+globalj])
loca_lmat[locali*n + globalj] = temp;
}
}
free(ro_wk);
}
int Owner(int k, int p, int n) { return k/(n/p);
}


void Copy_row(int loca_lmat[], int n, int p, int ro_wk[], int k) { int j;
int local_k = k % (n/p); for (j = 0; j < n; j++)
ro_wk[j] = loca_lmat[local_k*n + j];
}
