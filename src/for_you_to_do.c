#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 129;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, maxind;
    double max;


    for(i = 0; i < n - 1; i++){
        maxind = i;
        max = fabs(A[i*n + i]);
        int t;

        // pivoting the matrix to find the maximum number
        for(t = i+1; t < n; t++){
            if(fabs(A[t*n + i]) > max){
                maxind = t;
                max = fabs(A[t*n + i]);
            }
        }

        if(max == 0){
            printf("LU factoration failed: coefficient matrix is singular");
            return -1;
        }
        else{

            //case when need to swap the row
            if(maxind != i){
                //swap pivoting information
                int temps= ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = temps;

                //swap row for matrix method 1
                // int j;
                // for(j = 0; j < n; j++){
                //     double k;
                //     k = A[i * n + j];
                //     A[i * n + j] = A[maxind * n + j];
                //     A[maxind * n + j] = k;
                // }

                //swap row method 2 -- need to test which one is faster
                double trow[n];
                memcpy(trow, A + i * n, n*sizeof(double));
                memcpy(A + i * n, A + maxind * n, n*sizeof(double));
                memcpy(A + maxind * n, trow, n*sizeof(double));
            }

        }

        //update A(i+1:n, i+1:n)
        int j;
        for(j = i + 1; j <n;j++){
            A[j*n + i] = A[j*n + i] / A[i*n + i];
            int k;
            for(k =  i + 1; k < n; k++){
                A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
            }
        }
    }

    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i,j;
    double sum;
    /* add your code here */
    if(UPLO == 'L'){
        //method 1 local array
        double y[n];

        //method 2, alloc space
        // double *y = (double *)malloc(sizeof(double)*n);

        //compare: not much difference in time, so just use method 1
        
        y[0] = B[ipiv[0]];
        for(i = 1; i < n; i++ ){
           sum = 0;
           for(j = 0; j < i; j++ ){
               sum += y[j] * A[i*n + j];
           } 
           y[i] = B[ipiv[i]] - sum;
        }

        for(i = 0; i < n; i++){
            B[i] = y[i];
        }
        // free(y);
    }

    else if(UPLO == 'U'){
        // double x[n];
        double *x = (double *)malloc(sizeof(double)*n);
        int i,j;
        double sum;

        x[n-1] = B[n-1] / A[(n-1) * n + n - 1];

        for(i = n - 2; i >= 0; i-- ){
            sum = 0;
            for(j = i + 1; j < n; j++){
                sum += x[j] * A[i*n + j];
            }
            x[i] = (B[i] - sum) / A[i*n + i];
        }
        for(i = 0; i < n; i++){
            B[i] = x[i];
        }
        // free(x);
    }
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    int ic,jc,kc,i1,j1,k1;
    register int m;
    int block_size = 3;


    for(i1 = i; i1 < n; i1 += b){

        for(j1 = j;j1 < n; j1 += b){

            for(k1 = k; k1 < k + b; k1 += b){




                for (ic=i1; ic<(i1 + b); ic+=block_size){


                    for (jc=j1; jc<(j1 + b); jc+=block_size) {




                        //9 register used for matrix C
                        register double c00 = C[ic * n + jc];
                        register double c01 = C[ic * n + (jc + 1)];
                        register double c02 = C[ic * n + (jc + 2)];

                        register double c10 = C[(ic + 1) * n + jc];
                        register double c11 = C[(ic + 1) * n + (jc + 1)];
                        register double c12 = C[(ic + 1) * n + (jc + 2)];

                        register double c20 = C[(ic + 2) * n + jc];
                        register double c21 = C[(ic + 2) * n + (jc + 1)];
                        register double c22 = C[(ic + 2) * n + (jc + 2)];

                        // double test = c00 + c01 + c02 + c10 + c11 + c12 + c20 + c21 + c22;

                        // printf("test = %f", test);


                        //6 registers INIT for A and B matrix
                        register double a00;
                        register double a10;
                        register double a20;
                        // register double a30;

                        register double b00;
                        register double b01;
                        register double b02;
                        // register double b03;

                        for (kc=k1; kc<(k1 + b); kc+=block_size){


                            //use for debug

                            // if(kc >= n){
                            //     printf("kc error\n");
                            //     printf("i1 = %i, j1 = %i, k1 = %i, ic = %i, jc = %i, kc = %i\n",i1,j1,k1,ic,jc,kc);
                            //     return;
                            // }


                            for(m = 0; m < block_size; m++){


                                
                                a00 = A[ic * n + kc + m];
                                a10 = A[(ic + 1)*n + kc + m];
                                a20 = A[(ic + 2)*n + kc + m];

                                b00 = B[(kc + m) * n + (jc)];
                                b01 = B[(kc + m) * n + (jc + 1)];
                                b02 = B[(kc + m) * n + (jc + 2)];

                                //Start doing the computing process
                                c00 -= a00 * b00;
                                c01 -= a00 * b01;
                                c02 -= a00 * b02;
                                c10 -= a10 * b00;
                                c11 -= a10 * b01;
                                c12 -= a10 * b02;
                                c20 -= a20 * b00;
                                c21 -= a20 * b01;
                                c22 -= a20 * b02;
                            }
                            

                        }
                        //Write back the value to matrix C
                        C[ic * n + jc] = c00;
                        C[ic * n + (jc + 1)] = c01;
                        C[ic * n + (jc + 2)] = c02;

                        C[(ic + 1) * n + jc] = c10;
                        C[(ic + 1) * n + (jc + 1)] = c11;
                        C[(ic + 1) * n + (jc + 2)] = c12;

                        C[(ic + 2) * n + jc] = c20;
                        C[(ic + 2) * n + (jc + 1)] = c21;
                        C[(ic + 2) * n + (jc + 2)] = c22;

                    }
                }

            }
        }
    }
    return;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{

    /* add your code here */

    int i,j,k,ic,t, maxind;
    double max;

    for(ic = 0; ic <n - 1;ic +=b){
        // if(ic > n){
        //     printf("error in ic\n\n\n");
        // }

        for(i = ic; i < ic+b ; i++){
            // if(i > n){
            //     printf("error in i\n\n\n");
            // }

            // pivoting the matrix to find the maximum number
            maxind = i;
            max = fabs(A[i*n + i]);
            for(t = i+1; t < n; t++){

                if(fabs(A[t*n + i]) > max){
                    maxind = t;
                    max = fabs(A[t*n + i]);
                }
            }

            
            if(max == 0){
                printf("LU factoration failed: coefficient matrix is singular");
                return -1;
            }
            else{

                //The case that need to swap the row
                if(maxind != i){
                    //swap pivoting information
                    int temps= ipiv[i];
                    ipiv[i] = ipiv[maxind];
                    ipiv[maxind] = temps;

                    //swap row for matrix method 1
                    // int j;
                    // for(j = 0; j < n; j++){
                    //     double k;
                    //     k = A[i * n + j];
                    //     A[i * n + j] = A[maxind * n + j];
                    //     A[maxind * n + j] = k;
                    // }

                    //swap row method 2 -- seem like not too much difference than method one, but this look a bit cleaner
                    double trow[n];
                    memcpy(trow, A + i * n, n*sizeof(double));
                    memcpy(A + i * n, A + maxind * n, n*sizeof(double));
                    memcpy(A + maxind * n, trow, n*sizeof(double));
                }

            }

            //update the lower triangle of A(ic:end , ic:end) and A(end+1:n , ic:end)
            for(j = i + 1; j <n;j++){
                // if(j > n){
                //     printf("error in j\n\n\n");
                // }

                A[j*n + i] = A[j*n + i] / A[i*n + i];

                //block version
                for(k = i + 1; k < ic + b; k++){

                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }

                //naive version - to test the top part of the code work
                // for(k = i + 1; k < n; k++){
                //     A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
                // }
            }
        }


        //update A(ic:end, end+1:n), basically same method as before, use the value store in A(ic:n, ic:end)
        register double total;
        //end = ic + b
        for(i = ic; i < ic + b; i++){
            // if(i > n){
            //     printf("error in i\n\n\n");
            // }

            for(j= ic + b;j < n;j++){
                // if(j > n){
                //     printf("error in j\n\n\n");
                // }

                total = 0;
                for(k = ic; k < i; k++){
                    // if(k > n){
                    //     printf("error in k\n\n\n");
                    // }

                    //naive version, abandon
                    // A[i*n - j] -= A[i*n + k] * A[k*n + j];

                    //new version, reduce access element
                    total += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= total;
            }
        }

        // update A(end + 1: n , end + 1 : n)
        // end = ic + b
        mydgemm(A, A, A,n, ic + b, ic + b, ic, b);
    }

    

    return 0;
}

