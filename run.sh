#!/bin/bash

#Compile
gcc -O2 ./Matrix/matrix_gen.c -o ./Matrix/matrix_gen
gcc ./SM/ge_sm_v1.c -o ./SM/ge_sm_v1 -lpthread
gcc ./SM/ge_sm_v2.c -o ./SM/ge_sm_v2 -lpthread
mpicc ./MPI/ge_mpi_v1.c -o ./MPI/ge_mpi_v1
mpicc ./MPI/ge_mpi_v2.c -o ./MPI/ge_mpi_v2

#Initialization
printBool=""
matrixSize=100
if [ $# -gt 0 ]; then
	if [ "$1"x = "-px" ]; then
		printBool=$1
	else
		matrixSize=$1
	fi
fi
if [ $# -gt 1 ]; then
	if [ "$2"x = "-px" ]; then
		printBool=$2
	else
		matrixSize=$2
	fi
fi
echo 'Do Gussian Elimination to a '$matrixSize'*'$matrixSize' matrix:'
matrixName=$matrixSize".matrix"

#Generate matrix
./Matrix/matrix_gen $matrixSize > $matrixName


for parallel in 1 2 4 6 8 10
do
solution_sm=$matrixSize"."$parallel"x.sm"
solution_mpi=$matrixSize"."$parallel"x.mpi"

echo '*** With '$parallel'x parallel:'
t_start=$(($(date +%s%N)/1000000))
./SM/ge_sm_v1 $parallel $printBool < $matrixName > $solution_sm"_v1"
t_end=$(($(date +%s%N)/1000000))
t_diff=$(($t_end - $t_start))
echo 'Shared memory V1 takes '$t_diff' msec.'

t_start=$(($(date +%s%N)/1000000))
./SM/ge_sm_v2 $parallel $printBool < $matrixName > $solution_sm"_v2"
t_end=$(($(date +%s%N)/1000000))
t_diff=$(($t_end - $t_start))
echo 'Shared memory V2 takes '$t_diff' msec.'

t_start=$(($(date +%s%N)/1000000))
mpiexec -np $(($parallel+1)) ./MPI/ge_mpi_v1 $printBool < $matrixName > $solution_mpi"_v1"
t_end=$(($(date +%s%N)/1000000))
t_diff=$(($t_end - $t_start))
echo 'Message passing interface V1 takes '$t_diff' msec.'

t_start=$(($(date +%s%N)/1000000))
mpiexec -np $parallel ./MPI/ge_mpi_v2 $printBool < $matrixName > $solution_mpi"_v2"
t_end=$(($(date +%s%N)/1000000))
t_diff=$(($t_end - $t_start))
echo 'Message passing interface V2 takes '$t_diff' msec.'
echo ''

if [ "$printBool"x = "x" ]; then
	rm $solution_sm"_v1"
	rm $solution_sm"_v2"
	rm $solution_mpi"_v1"
	rm $solution_mpi"_v2"
fi
done

#remove
rm $matrixName
