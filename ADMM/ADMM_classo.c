

#include "math.h"
#include "mex.h"
#include "blas.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    /* Pass inputs from Matlab and define constants */
	double *Y, *X, *C, *C_t, *z, *w, *rho, *ep_r, *ep_a, *lambda; 
    // Y = X'*Y+rho*(z_o-u_o)+rho*C'*(-w_o-t_o); X = (X'*X+rho*eye(n)+rho*C'*C)^(-1) in the formula.
	double *maxit, *print;
    double one = 1.0;
    double negOne = -1.0;
    double zero = 0.0;
    
    int iter = 0;
       
    double ep_p = zero; 
    double ep_d = zero;
    double s = zero, r_1 = zero, r_2 = zero; // intermediate variables to compute r2, s2 
    double r2 = one, s2 = one;
    
    int i, j;
//     ptrdiff_t *ipiv;
    
    ptrdiff_t l, n, m; // recording matrix dimensions
    
    ptrdiff_t info; // Used in blas lapack functions
    char *ch_N = "N"; 
    char *ch_U = "U";
    
    // geting arguments (as pointers) using mxGetPr
    Y = mxGetPr(prhs[0]);
    X = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    C_t = mxGetPr(prhs[3]);
    z = mxGetPr(prhs[4]);
    w = mxGetPr(prhs[5]);
    rho = mxGetPr(prhs[6]);
    ep_r = mxGetPr(prhs[7]);
    ep_a = mxGetPr(prhs[8]);
    lambda = mxGetPr(prhs[9]);
//     u = mxGetPr(prhs[10]);
//     t = mxGetPr(prhs[11]);
    maxit = mxGetPr(prhs[10]);
    print = mxGetPr(prhs[11]);
    
    /* dimensions of in put matrices */
    l = mxGetN(prhs[2]);
    n = mxGetM(prhs[2]);
    m = mxGetN(prhs[0]);
    
    double *z_o = z;
    double *w_o = w;
    mxArray *u, *t;
    double *u_o, *t_o;
    u = mxCreateDoubleMatrix(l, 1, mxREAL);
    u_o = mxGetPr(u);
    t = mxCreateDoubleMatrix(n, 1, mxREAL);
    t_o = mxGetPr(t);
    
//     for(i=0; i<10; i++){
// //         printf("*** %f Y*", Y[i]); 
//         printf("*** %f z* \n", z_o[i]);
//         printf("*** %d \n",i);
// 
// //         printf("*** %f w*", w_o[i]);        
//     }
    
//     Awork = mxCreateDoubleMatrix(m, p, mxREAL);
//     A2 = mxGetPr(Awork);
//     plhs[0] = mxCreateDoubleMatrix(p, n, mxREAL);
//     B2 = mxGetPr(plhs[0]);
//     memcpy(A2, A, m*p*mxGetElementSize(prhs[0]));
//     memcpy(B2, B, p*n*mxGetElementSize(prhs[1]));
    
    // intermediate variables used in the algorithm
	double *temp1 = (double *) mxMalloc(l*sizeof(double));  // rho*(z_o-u_o)
	double *temp2 = (double *) mxMalloc(n*sizeof(double));  // rho*(-w_o-t_o)
	double *temp3 = (double *) mxMalloc(l*sizeof(double));  // rho*C'*(-w_o-t_o)
	double *temp4 = (double *) mxMalloc(l*sizeof(double));  // X'*Y+rho*(z_o-u_o)+rho*C'*(-w_o-t_o)
	double *temp5 = (double *) mxMalloc(n*sizeof(double));  // C*beta_now
	double *temp6 = (double *) mxMalloc(n*sizeof(double));  // w_n-w_o
	double *temp7 = (double *) mxMalloc(l*sizeof(double));  // C'*(w_n-w_o)
	double *beta = (double *) mxMalloc(l*sizeof(double));   
    
    double temp82_1, temp82_2, temp92_1, temp92_2, temp8, temp9, temp10, temp12, temp122;
    
	double *temp11 = (double *) mxMalloc(l*sizeof(double));
   
//     double *beta_now;
//     int *iter_now;
//     double *beta_now = (double *) mxMalloc(l*sizeof(double));
//     int *iter_now = (int *) mxMalloc(1*sizeof(int));
//     double *r2_now = (double *) mxMalloc(1*sizeof(double));
//     double *s2_now = (double *) mxMalloc(1*sizeof(double));
//     double *ep_p_now = (double *) mxMalloc(1*sizeof(double));
//     double *ep_d_now = (double *) mxMalloc(1*sizeof(double));
//     double *z_now = (double *) mxMalloc(l*sizeof(double));
//     double *w_now = (double *) mxMalloc(n*sizeof(double));
//  
    double *z_n = (double *) mxMalloc(l*sizeof(double));
	double *w_n = (double *) mxMalloc(n*sizeof(double));
	double *u_n = (double *) mxMalloc(l*sizeof(double));
	double *t_n = (double *) mxMalloc(n*sizeof(double));

    /* create output matrix */
    plhs[0] = mxCreateDoubleMatrix(l, 1, mxREAL);
    beta = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(l, 1, mxREAL);
    z_n = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL);
    w_n = mxGetPr(plhs[2]);
    
//     plhs[3] = mxCreateDoubleMatrix(l, 1, mxREAL);
//     z_o = mxGetPr(plhs[3]);
//     
//     plhs[4] = mxCreateDoubleMatrix(n, 1, mxREAL);
//     w_o = mxGetPr(plhs[4]);
    
//     plhs[3] = mxCreateDoubleMatrix(l, 1, mxREAL);
//     u_n = mxGetPr(plhs[3]);
//     
//     plhs[4] = mxCreateDoubleMatrix(n, 1, mxREAL);
//     t_n = mxGetPr(plhs[4]);
    
//     plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
//     r2 = mxGetPr(plhs[2]);
//     
//     plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
//     s2 = mxGetPr(plhs[3]);
//     
//     plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
//     ep_p = mxGetPr(plhs[4]);
//     
//     plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
//     ep_d = mxGetPr(plhs[5]);
//     
//     plhs[2] = mxCreateDoubleMatrix(l, 1, mxREAL);
//     z_n = mxGetPr(plhs[2]);
//     
//     plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);
//     w_n = mxGetPr(plhs[3]);
//     
//     plhs[4] = mxCreateDoubleMatrix(l, 1, mxREAL);
//     beta_now = mxGetPr(plhs[4]);
//     
//     plhs[9] = mxCreateDoubleMatrix(n, 1, mxREAL);
//     temp5 = mxGetPr(plhs[9]);
    
//     plhs[8] = mxCreateDoubleMatrix(l, 1, mxREAL);
//     temp4 = mxGetPr(plhs[8]);
//     
//     plhs[9] = mxCreateDoubleMatrix(l, 1, mxREAL);
//     temp3 = mxGetPr(plhs[9]);
    
//     printf(" %f", s2);
//     printf(" %f", r2);
//     printf(" %f", ep_p);
//     printf(" %f", ep_d);
//     printf(" %f", *maxit);
//     printf(" %f", C[1]);
//     printf(" %f", C[2]);
//     printf(" %f", C[3]);
//     printf(" %f", C[4]);
    
    while((r2>ep_p || s2>ep_d) && iter<*maxit){
        
        for(i=0; i<l; i++){
            temp1[i] = *rho*(z_o[i]-u_o[i]);
        }
        
//         for(i=0; i<10; i++){
//             printf("*** %f z_o*", z_o[i]);
// //             printf("*** %f u_o*", u_o[i]);
// //             printf("*** %f rho*", rho);
// //             printf("*** %f temp1*", temp1[i]);
//         }
        
        for(i=0; i<n; i++){
            temp2[i] = *rho*(-w_o[i]-t_o[i]);
        }
        
        dgemm_(ch_N, ch_N, &l, &m, &n, &one, C_t, &l, temp2, &n, &zero, temp3, &l);
        
        for(i=0; i<l; i++){
            temp4[i] = Y[i]+temp1[i]+temp3[i];
        }
        
//         for(i=0; i<10; i++){
//             printf("*** %f temp4*", temp1[i]+temp3[i]);
//         }

        dgemm_(ch_N, ch_N, &l, &m, &l, &one, X, &l, temp4, &l, &zero, beta, &l); // beta_now
        
//         for(i=0; i<l; i++){
//             beta_now[i] = beta[i];
//         }
        dgemm_(ch_N, ch_N, &n, &m, &l, &one, C, &n, beta, &l, &zero, temp5, &n); //temp5
        
        // update z, u
        for(i=0; i<l; i++){
            z_n[i] = (fabs(beta[i]+u_o[i]) - *lambda / *rho)*((fabs(beta[i]+u_o[i])- *lambda / *rho)>0)*(((beta[i]+u_o[i])>0)-((beta[i]+u_o[i])<0));
            u_n[i] = u_o[i]+beta[i]-z_n[i];
        }

        // update w, t
        for(i=0; i<n; i++){
            if(-temp5[i]-t_o[i]>0){
                w_n[i] = -temp5[i]-t_o[i];
            } else {
                w_n[i] = 0;
            }
            t_n[i] = t_o[i]+temp5[i]+w_n[i];
        }
        
//         for(i=0; i<l; i++){
//             u_n[i] = u_o[i]+beta[i]-z_n[i];
//         }
        
//         for(i=0; i<n; i++){
//             t_n[i] = t_o[i]+temp5[i]+w_n[i];
//         }
        
        // compute r2, s2, ep_p, ep_d
        r_1 = 0.0;
        for(i=0; i<l; i++){
            r_1 += pow((beta[i]-z_n[i]),2);
        }
        
        r_2 = 0.0;
        for(i=0; i<n; i++){
            r_2 += pow((temp5[i]+w_n[i]),2);
        }
        r2 = sqrt(r_1+r_2);
        
        for(i=0; i<n; i++){
            temp6[i] = w_n[i]-w_o[i];
        }
        
        dgemm_(ch_N, ch_N, &l, &m, &n, &one, C_t, &l, temp6, &n, &zero, temp7, &l);

        s = 0.0;
        for(i=0; i<l; i++){
            s += pow((z_n[i]-z_o[i]-temp7[i]), 2);
//             printf("%f||%f||%f||%f||\n", s, z_n[i],z_o[i],temp7[i]);
        }
        
//         for(i=1; i<10; i++){
//             printf("||%f||%f||%f||\n", z_n[i],z_o[i],temp7[i]);
//         }
//         printf("s= %f\n",s);
        s2 = (*rho)*sqrt(s);
//         printf("s2= %f\n",s2);
        
        temp82_1 = 0.0;
        temp92_1 = 0.0;
        for(i=0; i<l; i++){
            temp82_1 += pow(beta[i],2);
            temp92_1 += pow(z_n[i],2);
        }
        
        temp82_2 = 0.0;
        temp92_2 = 0.0;
        for(i=0; i<n; i++){
            temp82_2 += pow(temp5[i],2);
            temp92_2 += pow(w_n[i],2);
        }
        
        temp8 = sqrt(temp82_1+temp82_2);
        temp9 = sqrt(temp92_1+temp92_2);
        if(temp8>temp9){
            temp10 = temp8;
        } else {
            temp10 = temp9;
        }
        
        ep_p = sqrt(n+l)*(*ep_a)+(*ep_r)*temp10;
        dgemm_(ch_N, ch_N, &l, &m, &n, &one, C_t, &l, t_n, &n, &zero, temp11, &l);
        
        temp12 = 0.0;
        for(i=0; i<l; i++){
            temp12 += pow(u_n[i],2)+2*u_n[i]*temp11[i]+pow(temp11[i],2);
        }
        temp122 = sqrt(temp12);
        ep_d = sqrt(l)*(*ep_a) + (*ep_r)*(*rho)*temp122;
        
        // updating z_o, u_o, w_o, t_o;
        for(i=0; i<n; i++){
            w_o[i] = w_n[i];
            t_o[i] = t_n[i];
        }
        
        for(i=0; i<l; i++){
            z_o[i] = z_n[i];
            u_o[i] = u_n[i];
        }
//         z_o = z_n;
//         w_o = w_n;
//         u_o = u_n;
//         t_o = t_n;
        
//         SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
//         dsytrs_(ch_U, &l, &m, X, &l, ipiv, temp4, &l, &info);
//         DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
//         DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
//         dgesv_(&l, &m, X, &l, ipiv, temp4, &l, &info);
//         dposv_(ch_U, &l, &m, X, &l, temp4, &l, &info);
//         DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
//         printf("info %d \n", info);
//         
//         
//         dgemm_(ch_N, ch_N, &n, &m, &l, &one, C, &n, temp3, &n, &zero, temp4, &n);
       
//         z_n = (abs(beta_now+u_o)-lambda/rho).*((abs(beta_now+u_o)-lambda/rho)>0).*sign(beta_now+u_o);
// (fabs(temp3[i]+u_o[i]) - *lambda / *rho)*((fabs(temp3[i]+u_o[i])- *lambda / *rho)>0)*(((temp3[i]+u_o[i])>0)-((temp3[i]+u_o[i])<0))
        
//         dgemm_(chn, chn, &m, &n, &p, &one, A, &m, B, &p, &zero, C, &m);
//     DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
//         for(i=0; i<n; i++){
//             B[i] = Y[i]+*rho*(z[i]-u[i]);
//         }
        iter++;
//         printf(" %f", r2);
//         printf(" %f", s2);
//         printf(" %f", ep_p);
//         printf(" %f \n", ep_d);
    }
    if(*print==1){
        printf(" %d \n", iter);
    }
//        beta_now = (X'*X+rho*eye(n)+rho*C'*C)\(X'*Y+rho*(z_o-u_o)+rho*C'*(-w_o-t_o));
    
    /* Pass arguments to Fortran by reference */
    //dgemm(ch_N, ch_N, &l, &one, &n, &one, C, &l, temp2, &n, &zero, temp3, &l);
    
//     beta_now = beta;
//     *iter_now = iter;
//     r2_now = &r2;
//     s2_now = &s2;
//     ep_p_now = &ep_p;
//     ep_d_now = &ep_d;
//     z_now = z_n;
//     w_now = w_n;
    
//     printf("   %d  \n",l);
//     printf("   %d  \n",n);
//     printf("   %d  \n",m);
//     printf("   %f  \n",one);
//     printf("   %f  \n",zero);
    mxFree(temp1);
    mxFree(temp2);
    mxFree(temp3);
    mxFree(temp4);
    mxFree(temp5);
    mxFree(temp6);
    mxFree(temp7);
    mxFree(temp11);
}