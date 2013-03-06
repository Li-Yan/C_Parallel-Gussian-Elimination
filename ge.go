package main

import (
	"fmt"
	"math/rand"
	"time"
	"runtime"
)

func subtractSlice (ar1, ar2 []float32, fac float32, ch chan int){
	for i := 0; i < len(ar1); i++{
		ar2[i] -= ar1[i]*fac
	}
	ch <- 1
}

func ge_serial (M [][]float32){
	N := len(M)
	for p := 0; p < N; p++ {
		for r := p + 1; r < N; r++ {
			ratio := M[r][p] / M[p][p]
			M[r][p] = 0.0
			for c := p + 1; c < N; c++ {
				M[r][c] -= ratio * M[p][c]
			}
		}
	}	
}

func ge_parallel (M [][]float32){
	N := len(M)
	for p := 0; p < N; p++ {
		ch := make(chan int, N - p -1)
		for r := p + 1; r < N; r++ {
			ratio := M[r][p] / M[p][p]
			M[r][p] = 0.0
			go subtractSlice(M[p][p+1:N], M[r][p+1:N], ratio, ch)
		}
		for r := p + 1; r < N; r++ {
			<- ch
		}
	}
}

func initialze() {
	
}

func main() {
	N := 2000
	runtime.GOMAXPROCS(runtime.NumCPU())
	// allocate memory
	M := make([][]float32, N)
	for i := range M {
		M[i] = make([]float32, N)
	}

	var i, j int
	var diag float32
	diag = 2.0
	for i = 0; i < N; i++ {
		for j = 0; j < N; j++ {
			M[i][j] = rand.Float32()
			if i == j {
				M[i][j] += diag
			}
		}
	}

	t0 := time.Now()
	ge_serial(M)
	t1 := time.Now()

	fmt.Printf("Serial gaussian elimination took %v to run.\n", t1.Sub(t0))

	MM := make([][]float32, N)
	for i := range M {
		MM[i] = make([]float32, N)
	}

	for i = 0; i < N; i++ {
		for j = 0; j < N; j++ {
			MM[i][j] = rand.Float32()
			if i == j {
				MM[i][j] += diag
			}
		}
	}

	t0 = time.Now()
	ge_parallel(MM)
	t1 = time.Now()

	fmt.Printf("Parallel gaussian elimination took %v to run.\n", t1.Sub(t0))

	// for i = 0; i < N; i++{
 //    	for j = 0; j < N; j++{
 //    		fmt.Printf("%f\t", MM[i][j])
 //    	}
 //    	fmt.Printf("\n")
 //    }
	
	M = nil
	MM = nil
}
