#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Define here your functions
void trmv_row(double *x, double *y, double **A, int n, int nthreads, omp_sched_t opt)
{
    int i,j;
    for( i=0; i<n; i++) y[i] = 0;
    omp_set_schedule(opt,0);
#pragma omp parallel for num_threads(nthreads) private(j) firstprivate(A,x,y)
    for( i=0; i<n; i++)
        for( j=0; j<i+1; j++)
            y[i] = y[i] + A[i][j] * x[j];
}
void trmv_col(double *x, double *y, double **A, int n, int nthreads, omp_sched_t opt)
{
    int i,j;
    double *y_tmp;
    y_tmp = (double *) malloc(sizeof(double)*n);
    for( i=0; i<n; i++) y_tmp[i] = 0;
    omp_set_schedule(opt,0); 
#pragma omp parallel for num_threads(nthreads) private(i) firstprivate(A,x,y_tmp)
    for( j=0; j<n; j++)
        for( i=j; i<n; i++)
            y_tmp[i] = y_tmp[i] + A[i][j] * x[j];

    //for( i=0; i<n; i++)
    //y[i]=y_tmp[n-i-1];
}
void best_seq(double *x, double *y, double *y_row, double *y_col, double **A, int n)
{
    int i;
    trmv_col(x, y_col, A, n, 1,1);
    trmv_row(x, y_row, A, n, 1,1);
    for( i=0; i<n; i++) y[i] = y_row[i]+y_col[i];
}
void best_par(double *x, double *y, double *y_row, double *y_col, double **A, int n, int nthreads, omp_sched_t opt)
{
    int i;
    trmv_col(x, y_col, A, n, nthreads,opt);
    trmv_row(x, y_row, A, n, nthreads,opt);
    omp_set_schedule(opt,0); // set cheunk to 0 => default
#pragma omp parallel for num_threads(nthreads)
    for( i=0; i<n; i++) y[i] = y_row[i]+y_col[i];
}
int main( void )
{
    int i,j,count=0;
    //omp_sched_t sch_static = 1;
    //omp_sched_t sch_dynamic = 2;
    //omp_sched_t opt;

    int nsizes = 10, 
        sizes[nsizes],
        n, rep, nreps = 100; 
        //sizes[nsizes] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000},
        for(i=0;i<nsizes;i++) sizes[i] = (i+1)*100;
        double t0, t1, tt, tt_seq;
        int nthreads;

        // Create/Allocate and initialize A, x, and y here
        double *x, *y, *y_row, *y_col, **A;
        x = (double *) malloc(sizeof(double)*sizes[nsizes-1]);
        y = (double *) malloc(sizeof(double)*sizes[nsizes-1]);
        y_row = (double *) malloc(sizeof(double)*sizes[nsizes-1]);
        y_col = (double *) malloc(sizeof(double)*sizes[nsizes-1]);
        A = malloc(sizeof *A * sizes[nsizes-1]);
        if (A) for (i = 0; i < sizes[nsizes-1]; i++) A[i] = malloc(sizeof *A[i] * sizes[nsizes-1]);
        for( i=0; i<sizes[nsizes-1]; i++) 
            for( j=0; j<i+1; j++) // only lower triangle of A
            {
                A[i][j] = i*j+i-j;
                x[i] = 0.0;
                y_row[i] = 0.0;
                y_col[i] = 0.0;
            }

        printf("Version                        | #size | #threads | #Time (s) | Speedup \n"); // Speedup wrt the base case
        printf("------------------------------------------------------------- \n");

        for (i = 0; i < nsizes; i++)
        {
            n = sizes[i];
            /*
               t0 = omp_get_wtime();
            // Here call your base case
            nthreads = 1;
            for (rep = 0; rep < nreps; rep++)
            best_seq(x,y,y_row,y_col,A,n);
            t1 = omp_get_wtime();
            tt_seq = t1 - t0;
            printf("Sequential                     | %5d | %8d | %9.5f | %7.2f \n", n, nthreads, tt_seq, tt_seq/tt_seq);
            */

            // Add your code here! (based on the code in lines 27-33 above)
            // for each function
            //    for each schedule
            //       for each number of threads in [1, 2, 4, 8]
            // Total of 160 function calls and printf statements
            //
            // Use the appropriate OpenMP functions to set number of threads and schedules
            //   to reduce code replication!!!
            // Here call your base case
            /* -----------  1 thread  ---------- */

            t0 = omp_get_wtime();
            nthreads = 1;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt_seq = t1 - t0;
            printf("\033[22;35mSequential   row (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt_seq);

            t0 = omp_get_wtime();
            nthreads = 1;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt_seq = t1 - t0;
            printf("\033[22;35mSequential   row (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt_seq);

            t0 = omp_get_wtime();
            nthreads = 1;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt_seq = t1 - t0;
            printf("\033[22;35mSequential   col (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt_seq);

            t0 = omp_get_wtime();
            nthreads = 1;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt_seq = t1 - t0;
            printf("\033[22;35mSequential   col (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt_seq);
            t0 = omp_get_wtime();


            /* -----------  2 thread  ---------- */
            nthreads = 2;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;33mParallelized row (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 2;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;33mParallelized row (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 2;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;33mParallelized col (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);
            
            t0 = omp_get_wtime();
            nthreads = 2;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;33mParallelized col (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            /* -----------  4 thread  ---------- */
            t0 = omp_get_wtime();
            nthreads = 4;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;32mParallelized row (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 4;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;32mParallelized row (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);
            
            t0 = omp_get_wtime();
            nthreads = 4;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;32mParallelized col (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 4;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;32mParallelized col (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            /* -----------  8 threads  ---------- */
            t0 = omp_get_wtime();
            nthreads = 8;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;31mParallelized row (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 8;
            for (rep = 0; rep < nreps; rep++) trmv_row(x, y_row, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;31mParallelized row (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 8;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,1);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_static);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;31mParallelized col (static)          | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);

            t0 = omp_get_wtime();
            nthreads = 8;
            for (rep = 0; rep < nreps; rep++) trmv_col(x, y_col, A, n, nthreads,2);
                //best_par(x,y,y_row,y_col,A,n,nthreads,sch_dynamic);
            t1 = omp_get_wtime(); tt = t1 - t0;
            printf("\033[22;31mParallelized col (dynamic)         | %5d | %8d | %9.5f | %7.2f \033[0m \n", n, nthreads, tt, tt_seq/tt);
            count = count + 16;
        }

        printf("------------------------------------------------------------- \n");
        //    for( i=0; i<n; i++) printf("y[%d]=%f\n",i,y[i]);
        printf("Print counter: %d\n", count);

        free(x);
        free(y);
        free(y_row);
        free(y_col);
        free(A);
        return 0;
}

