#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

double **matrix;	//shared memory
int matrixSize;
int num_threads;
pthread_barrier_t barr, all_done;

inline int Max_Int(int x, int y){
	return x ^ ((x ^ y) & -(x < y));
}

double** Read_Matrix(int *MatrixSize)
{
	int i, j;
	fscanf(stdin, "%d", MatrixSize);

	matrix = (double **) malloc(*MatrixSize * sizeof(double*));
	for (i = 0; i < *MatrixSize; i++)
	{
		matrix[i] = (double *) malloc(*MatrixSize * sizeof(double));
		for (j = 0; j < *MatrixSize; j++)
		{
			fscanf(stdin, "%lf", &matrix[i][j]);
		}
	}
	return matrix;
}

void Print_Matrix(const int MatrixSize)
{
	int i, j;
	for (i = 0; i < MatrixSize; i++)
	{
		for (j = 0; j < MatrixSize - 1; j++)
		{
			printf("%.3f ", matrix[i][j]);
		}
		printf("%.3f\n", matrix[i][MatrixSize - 1]);
	}
}

void Free_Matrix(const int MatrixSize)
{
	int i;
	for (i = 0; i < MatrixSize; ++i)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void* Gaussian_Elimination_SM_Thread(void* tn)
{
	long tid = (long)tn;
	int p;
	int block_size = matrixSize / num_threads;
	int rem = matrixSize % num_threads;
	int start_row = tid * block_size;
	int end_row = (tid + 1) * block_size - 1;
	if (tid == (num_threads-1))
	{
		end_row += rem;
	}	
	int r, c;
	double ratio;
	for (p = 0; p < matrixSize; p++)
	{
		for (r = Max_Int(p + 1, start_row); r <= end_row; r++)
		{
			ratio = matrix[r][p] / matrix[p][p];
			for (c = p; c < matrixSize; c++){
				matrix[r][c] -= ratio * matrix[p][c];
			}
		}
		pthread_barrier_wait(&barr);
	}
	pthread_barrier_wait(&all_done);
	pthread_exit(NULL);
}

void Gaussian_Elimination_SM()
{
	int j, k;
	long t;
	double ratio;
	pthread_t threads[num_threads];
	pthread_barrier_init(&barr, NULL, num_threads);
	pthread_barrier_init(&all_done, NULL, num_threads + 1);

	for (t = 0; t < num_threads; t++)
	{
		pthread_create(&threads[t],
						NULL, 
						Gaussian_Elimination_SM_Thread, 
						(void *) t);
	}
	pthread_barrier_wait(&all_done);
	return;
}

int main(int argc, char const *argv[])
{
	short print_bool = 0;

	if (argc > 2)
	{
		print_bool = strcmp(argv[2], "-p") == 0;
	}

	num_threads = atoi(argv[1]);
	Read_Matrix(&matrixSize);

	if (print_bool)
	{
		printf("******************** Original Matrix ********************\n");
		Print_Matrix(matrixSize);
	}

	Gaussian_Elimination_SM();

	if (print_bool)
	{
		printf("\n******************** New Matrix ********************\n");
		Print_Matrix(matrixSize);
	}

	Free_Matrix(matrixSize);
	return 0;
}
