#include <stdio.h>
#include <stdlib.h>

#define DIAG 2.0

// create the array
// the entries are DIAG*1(i == j) + rand(0,1)
double** Init_Matrix(int MatrixSize)
{
	int i, j;
	double ** matrix;
	matrix = (double**) malloc(MatrixSize * sizeof(double*));
	for (i = 0; i < MatrixSize; i++) {
		matrix[i] = (double*) malloc(MatrixSize * sizeof(double));
		for (j = 0; j < MatrixSize; j++){
			matrix[i][j] = ((double) rand()) / ((double) RAND_MAX);
			if (i == j){
				matrix[i][j] += DIAG;
			}
		}
	}
	return matrix;
}

void Print_Matrix(double **Matrix, const int MatrixSize)
{
	int i, j;
	for (i = 0; i < MatrixSize; i++)
	{
		for (j = 0; j < MatrixSize - 1; j++)
		{
			printf("%f ", Matrix[i][j]);
		}
		printf("%f\n", Matrix[i][MatrixSize - 1]);
	}
}

void Free_Matrix(double **Matrix, const int MatrixSize)
{
	int i;
	for (i = 0; i < MatrixSize; ++i)
	{
		free(Matrix[i]);
	}
	free(Matrix);
}

int main(int argc, char const *argv[])
{
	double **matrix;
	int matrixSize = atoi(argv[1]);

	printf("%d\n", matrixSize);
	matrix = Init_Matrix(matrixSize);
	Print_Matrix(matrix, matrixSize);
	Free_Matrix(matrix, matrixSize);
	return 0;
}