#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MANAGER_PROC			0
#define MATRIXSIZE_TAG 		227800000
#define OFFSET_MATRIX_TAG	700000000

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

double** Gaussian_Elimination_Sub(double *BaseRowMatrix, double **SubMatrix, const int RowNum, 
	const int MatrixSize, const int BaseRow)
{
	int i, j;
	double ratio;

	for (i = 0; i < RowNum; i++)
	{
		ratio = SubMatrix[i][BaseRow] / BaseRowMatrix[BaseRow];
		for (j = BaseRow; j < MatrixSize; j++)
		{
			SubMatrix[i][j] -= ratio * BaseRowMatrix[j];
		}
	}
	return SubMatrix;
}

int main(int argc, char **argv)
{
	int ierr, num_procs, num_worker, my_id, recv_id;
	int i, j, base_row, remain_row, start_row, end_row, block_size;
	int size_tag, base_row_tag, sub_matrix_tag, finish_tag, num_tags = 4;
	MPI_Status status;

	double **matrix, *matrix_base_row, **matrix_sub, dd;	//unshared
	int matrixSize, matrix_subSize;

	short print_bool = 0;

	if (argc > 1)
	{
		print_bool = strcmp(argv[1], "-p") == 0;
	}

	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	num_worker = num_procs - 1;

	if (my_id == 0)
		//Main Process
	{
		matrix = Read_Matrix(&matrixSize);
		if (print_bool)
			{
			printf("******************** Original Matrix ********************\n");
			Print_Matrix(matrix, matrixSize);
			printf("\n");
		}

		for (recv_id = 1; recv_id < num_procs; recv_id++)
		{
			ierr = MPI_Send(&matrixSize, 1 , MPI_INT, recv_id, MATRIXSIZE_TAG, MPI_COMM_WORLD);
		}

		for (base_row = 0; base_row < matrixSize - 1; base_row++)
		{
			//Tag initialization
			size_tag = base_row * num_tags + 0;
			base_row_tag = base_row * num_tags + 1;
			sub_matrix_tag = base_row * num_tags + 2;
			finish_tag = base_row * num_tags + 3;

			//Distribution initialization
			remain_row = matrixSize - base_row - 1;
			block_size = remain_row / num_worker;
			if (remain_row % num_worker != 0) block_size++;
			end_row = base_row;

			for (recv_id = 1; recv_id < num_procs; recv_id++)
				//Give works to other processes to do
			{
				start_row = end_row + 1;
				if (start_row < matrixSize)
					//There is work for this process
				{
					end_row = start_row + block_size - 1;
					if (end_row > matrixSize - 1) 
					{
						end_row = matrixSize - 1;
					}
					matrix_subSize = end_row - start_row + 1;
					ierr = MPI_Send(&matrix_subSize, 1 , MPI_INT, recv_id, size_tag, 
						MPI_COMM_WORLD);

					//Send base row & sub_matrix
					ierr = MPI_Send(&(matrix[base_row][0]), 2 * matrixSize, MPI_INT, recv_id, 
						base_row_tag, MPI_COMM_WORLD);
					for (i = 0; i < matrix_subSize; ++i)
					{
						ierr = MPI_Send(&(matrix[i + start_row][0]), 2 * matrixSize, MPI_INT, recv_id, 
							sub_matrix_tag + i + OFFSET_MATRIX_TAG, MPI_COMM_WORLD);
					}

					//Receive updated sub_matrix
					for (i = 0; i < matrix_subSize; ++i)
					{
						ierr = MPI_Recv(&(matrix[i + start_row][0]), 2 * matrixSize, MPI_INT, recv_id, 
							sub_matrix_tag + i + OFFSET_MATRIX_TAG, MPI_COMM_WORLD, &status);
					}
				}
				else
					//There is no work for him
				{
					matrix_subSize = 0;
					ierr = MPI_Send(&matrix_subSize, 1 , MPI_INT, recv_id, size_tag, 
						MPI_COMM_WORLD);
				}
			}

			for (recv_id = 1; recv_id < num_procs; recv_id++)
				//Wait for each process finish this turn
			{
				ierr = MPI_Recv(NULL, 0, MPI_INT, recv_id, finish_tag, 
					MPI_COMM_WORLD, &status);
			}
		}

		if (print_bool)
		{
			printf("******************** New Matrix ********************\n");
			Print_Matrix(matrix, matrixSize);
		}
		Free_Matrix(matrix, matrixSize);
	} 
	else 
		//Child Process
	{
		ierr = MPI_Recv(&matrixSize, 1, MPI_INT, MANAGER_PROC, MATRIXSIZE_TAG, 
			MPI_COMM_WORLD, &status);

		for (base_row = 0; base_row < matrixSize - 1; base_row++)
		{
			//Tag initialization
			size_tag = base_row * num_tags + 0;
			base_row_tag = base_row * num_tags + 1;
			sub_matrix_tag = base_row * num_tags + 2;
			finish_tag = base_row * num_tags + 3;

			remain_row = matrixSize - base_row - 1;

			ierr = MPI_Recv(&matrix_subSize, 1, MPI_INT, MANAGER_PROC, size_tag, 
				MPI_COMM_WORLD, &status);
			if (matrix_subSize > 0)
			{
				//Receive base row
				matrix_base_row = (double *) malloc((matrixSize) * sizeof(double));
				ierr = MPI_Recv(&(matrix_base_row[0]), 2 * matrixSize, MPI_INT, MANAGER_PROC, 
					base_row_tag, MPI_COMM_WORLD, &status);
				
				//Receive sub_matrix
				matrix_sub = (double **) malloc(matrix_subSize * sizeof(double *));
				for (i = 0; i < matrix_subSize; ++i)
				{
					matrix_sub[i] = (double *) malloc(matrixSize * sizeof(double));
					ierr = MPI_Recv(&(matrix_sub[i][0]), 2 * matrixSize, MPI_INT, MANAGER_PROC, 
						sub_matrix_tag + i + OFFSET_MATRIX_TAG, MPI_COMM_WORLD, &status);
				}

				//Elimination
				matrix_sub = Gaussian_Elimination_Sub(matrix_base_row, matrix_sub, 
					matrix_subSize, matrixSize, base_row);

				//Send updated sub_matrix
				for (i = 0; i < matrix_subSize; ++i)
				{
					ierr = MPI_Send(&(matrix_sub[i][0]), 2 * matrixSize, MPI_INT, MANAGER_PROC, 
						sub_matrix_tag + i + OFFSET_MATRIX_TAG, MPI_COMM_WORLD);
				}

				//Free temporary memory
				free(matrix_base_row);
				for (i = 0; i < matrix_subSize; ++i)
				{
					free(matrix_sub[i]);
				}
				free(matrix_sub);
			}

			ierr = MPI_Send(NULL, 0 , MPI_INT, MANAGER_PROC, finish_tag, 
				MPI_COMM_WORLD);
		}
	}

  ierr = MPI_Finalize();

  return 0;
}
