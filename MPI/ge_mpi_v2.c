#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MANAGER_PROC			0
#define MATRIXSIZE_TAG 		227800000
#define FINISH_TAG 				600000000

double** Read_Matrix(int *MatrixSize)
{
	int i, j;
	double **matrix;
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

void Print_Matrix(double **Matrix, const int MatrixSize)
{
	int i, j;
	for (i = 0; i < MatrixSize; i++)
	{
		for (j = 0; j < MatrixSize - 1; j++)
		{
			printf("%.3f ", Matrix[i][j]);
		}
		printf("%.3f\n", Matrix[i][MatrixSize - 1]);
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

double** Gaussian_Elimination_Sub(double *Base_Row, double **Sub_Matrix, const int MatrixSize,
	const int Current_Row_Point, const int Start_Row_Point, const int End_Row_Point)
{
	int i, j;
	double ratio;

	for (i = Start_Row_Point; i <= End_Row_Point; i++)
	{
		ratio = Sub_Matrix[i][Current_Row_Point] / Base_Row[Current_Row_Point];
		for (j = Current_Row_Point; j < MatrixSize; j++)
		{
			Sub_Matrix[i][j] -= ratio * Base_Row[j];
		}
	}
	return Sub_Matrix;
}

int* Cal_All_Start_Row(const int MatrixSize, const int Num_Procs)
{
	int *start_row = (int *) malloc((Num_Procs + 1) * sizeof(int));
	int end_row;
	int block_size = MatrixSize / Num_Procs;
	int remaining_row_num = MatrixSize % Num_Procs;
	int i;

	end_row = -1;
	for (i = 0; i < Num_Procs; i++)
	{
		start_row[i] = end_row + 1;
		end_row = start_row[i] + block_size - 1;
		if (remaining_row_num > 0)
		{
			end_row++;
			remaining_row_num--;
		}
	}
	start_row[Num_Procs] = MatrixSize;
	return start_row;

}

int Find_Current_Row_Sender(int *Start_Row, const int Current_Row_Point, const int Num_Procs)
{
	int i;
	for (i = 1; i <= Num_Procs; i++)
	{
		if (Start_Row[i] > Current_Row_Point) break;
	}
	return i - 1;
}

int main(int argc, char **argv)
{
	int ierr, num_procs, num_worker, my_id, recv_id, send_id;
	int i, j, m, *start_row, row_num, block_size;
	int current_row_point;
	int size_tag, base_row_tag, sub_matrix_tag, finish_tag, num_tags = 4;
	MPI_Status status;

	double **sub_matrix, *base_row, dd;	//unshared
	int matrixSize, matrix_subSize;

	short print_bool = 0;

	if (argc > 1)
	{
		print_bool = strcmp(argv[1], "-p") == 0;
	}

	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	//First process read the matrix and braodcase the matrix size
	if (my_id == MANAGER_PROC)
	{
		sub_matrix = Read_Matrix(&matrixSize);
		if (print_bool)
		{
			printf("******************** Original Matrix ********************\n");
			Print_Matrix(sub_matrix, matrixSize);
			printf("\n");
		}
		
		//Broadcast matrix size
		for (recv_id = 1; recv_id < num_procs; recv_id++)
		{
			ierr = MPI_Send(&matrixSize, 1 , MPI_INT, recv_id, MATRIXSIZE_TAG, MPI_COMM_WORLD);
		}
	} 
	else
	{
		ierr = MPI_Recv(&matrixSize, 1, MPI_INT, MANAGER_PROC, MATRIXSIZE_TAG, MPI_COMM_WORLD, &status);
	}

	start_row = Cal_All_Start_Row(matrixSize, num_procs);

	//If the process number is greater than matrix size, some process has not job
	if (start_row[my_id] >= matrixSize)
	{
		free(start_row);
		ierr = MPI_Finalize();
		return 0;
	}
	for (i = num_procs - 1; i >= 0; i--)
	{
		if (start_row[i] < matrixSize) break;
	}
	num_procs = i + 1;
	row_num = start_row[my_id + 1] - start_row[my_id];

	//First process send matrix part to other processes
	if (my_id == MANAGER_PROC)
	{
		for (recv_id = 1; recv_id < num_procs; recv_id++)
		{
			for (i = start_row[recv_id]; i < start_row[recv_id + 1]; i++)
			{
				ierr = MPI_Send(&(sub_matrix[i][0]), matrixSize , MPI_DOUBLE, recv_id, 
					i - start_row[recv_id], MPI_COMM_WORLD);
			}
		}

		//Synchronize submatrix transmission
		for (recv_id = 1; recv_id < num_procs; recv_id++)
		{
			ierr = MPI_Recv(NULL, 0, MPI_INT, recv_id, FINISH_TAG, MPI_COMM_WORLD, &status);
		}
	}
	else 
	{
		sub_matrix = (double **) malloc(row_num * sizeof(double *));
		for (i = 0; i < row_num; i++)
		{
			sub_matrix[i] = (double *) malloc(matrixSize * sizeof(double));
			ierr = MPI_Recv(&(sub_matrix[i][0]), matrixSize , MPI_DOUBLE, MANAGER_PROC, 
				i, MPI_COMM_WORLD, &status);
		}
		ierr = MPI_Send(NULL, 0 , MPI_DOUBLE, MANAGER_PROC, FINISH_TAG, MPI_COMM_WORLD);
	}

	//Start Guass Elimination
	base_row = (double *) malloc(matrixSize * sizeof(double));
	for (current_row_point = 0; current_row_point < matrixSize - 1; current_row_point++)
	{
		if (start_row[my_id + 1] <= current_row_point)
			//This process has finished all its job
		{
			continue;
		}

		if ((start_row[my_id] <= current_row_point) && (start_row[my_id + 1] > current_row_point))
			//This process should broadcast the current row
		{
			for (recv_id = my_id + 1; recv_id < num_procs; recv_id++)
			{
				ierr = MPI_Send(&(sub_matrix[current_row_point - start_row[my_id]][0]), matrixSize , 
					MPI_DOUBLE, recv_id, current_row_point, MPI_COMM_WORLD);
			}

			if (current_row_point < start_row[my_id + 1] - 1)
			{
				Gaussian_Elimination_Sub(sub_matrix[current_row_point - start_row[my_id]], sub_matrix, 
					matrixSize, current_row_point, current_row_point - start_row[my_id] + 1, row_num - 1);
			}
		}
		else
		{
			send_id = Find_Current_Row_Sender(start_row, current_row_point, num_procs);
			ierr = MPI_Recv(&(base_row[0]), matrixSize , MPI_DOUBLE, send_id, 
				current_row_point, MPI_COMM_WORLD, &status);

			Gaussian_Elimination_Sub(base_row, sub_matrix, matrixSize, current_row_point, 0, row_num - 1);
		}
	}

	//Print result
	if (print_bool)
	{
		if (my_id != MANAGER_PROC)
		{
			ierr = MPI_Recv(NULL, 0, MPI_INT, my_id - 1, FINISH_TAG, MPI_COMM_WORLD, &status);
		}
		else
		{
			printf("******************** New Matrix ********************\n");
		}
		for (i = 0; i < row_num; i++)
		{
			for (j = 0; j < matrixSize - 1; j++)
			{
				printf("%.3f ", sub_matrix[i][j]);
			}
			printf("%.3f\n", sub_matrix[i][matrixSize - 1]);
		}
		if (my_id != num_procs - 1)
		{
			ierr = MPI_Send(NULL, 0, MPI_INT, my_id + 1, FINISH_TAG, MPI_COMM_WORLD);
		}
	}

	free(base_row);
	if (my_id = MANAGER_PROC) row_num = matrixSize;
	for (i = 0; i < row_num; i++)
	{
		free(sub_matrix[i]);
	}
	free(sub_matrix);
  ierr = MPI_Finalize();

  return 0;
}
