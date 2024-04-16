#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main ( int argc, char *argv[] )  
{
    int i, n;
    long long int *fib;
    int tid, nthreads, chunk;
    double start, end;
    const double phi = (1+sqrt(5))/2;
    if (argc != 2)
    {
        fprintf(stderr, "Usage %s <n>\n", argv[0]);
        exit(-1);
    }
    n = atoi(argv[1]);
    fib = malloc(n * sizeof(long long int));
    omp_set_num_threads(10); // run-time no. of threads spawn
    start = omp_get_wtime();
#pragma omp parallel private(tid) shared(nthreads)
{
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chunk = n/nthreads+1; // evenly distribute, +1 to compensate round-down
    // Compute
#pragma omp for schedule(static,chunk)
    for( i = 0; i < n; i++ )
    {
        //printf("[%d]:I am %d\n", i, tid); // check worksharing
        if(i<chunk*tid+2) // only use Binet's Formular at the fist two terms.
            fib[i] = (pow(phi,i+1)-(pow(-1*phi,-1*(i+1))))/sqrt(5); // Seeding
        else 
            fib[i] = fib[i-1] + fib[i-2]; // loop-carried dependencies
    //printf("Fib(%d): %lld\n", i, fib[i]);
    }
}
    end = omp_get_wtime();
    // Print out the n-th Fibonacci number
    printf("Fib(%d): %lld\n", n, fib[n-1]);
    printf("%d threads: %10.10f seconds\n", nthreads, end-start);

    free(fib);
    return 0;
}
