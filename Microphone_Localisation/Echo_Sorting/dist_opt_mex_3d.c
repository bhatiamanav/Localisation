/* [dist, D_estim, XY_coor_estim] = dist_opt_mex_3d(D, dim)
 * this program finds the optimal frobinius distance of a matrix to an EDM. It uses
 * the s-stress cost function and optimizes it with a coordinate optimization method.
 *
 * INPUTS:
 *       D                  : Input matrix for which we want to find the distance to an EDM.
 *                            Elements of D are squared distances.
 *       dim                : Ambient dimension.
 *
 * OUTPUTS:
 *       dist               : Frobenius distance between the input (D) and the 
 *                            closest EDM (D_estim).
 *       D_estim            : Optimal EDM closest to D with embedding dim.
 *       XY_coor_estim      : Estimated positions that yield D_estim.
 */

#include <math.h>
#include "mex.h"
#define pi 3.1415926535897932384626433


/* the function cubicroots finds the roots of a cubic polynomial with
 * coefficients a, b, c, d and puts the result in vector roots. */
static void cubicroots(double a, double	b, double c, double d, double *roots)
{
    // calculate p and q
    double p  = c/a - b*b/a/a/3. ;
    double q  = (2.*b*b*b/a/a/a - 9.*b*c/a/a + 27.*d/a) / 27. ;
    
    // compute the discriminant
    double DD = p*p*p/27. + q*q/4. ;
    
    double temp1;
    double temp2;
    double y1;
    double y2;
    double y3;
    double u;
    double v;
    double phi;
    double y2i;
    double y2r;
    
    if (DD < 0.)
    {
        // 3 real unequal roots -- use the trigonometric formulation
        phi = acos(-q/2./sqrt(fabs(p*p*p)/27.));
        temp1 = 2.*sqrt(fabs(p)/3.);
        y1 =  temp1*cos(phi/3.);
        y2 = -temp1*cos((phi+pi)/3.);
        y3 = -temp1*cos((phi-pi)/3.);
    }
    
    else
    {
        // 1 real root & 2 conjugate complex roots OR 3 real roots (not necessarily distinct)
        temp1 = -q/2. + sqrt(DD);
        temp2 = -q/2. - sqrt(DD);
        u = pow(fabs(temp1),(1./3.));
        v = pow(fabs(temp2),(1./3.));
        if (temp1 < 0.) u=-u;
        if (temp2 < 0.) v=-v;
        y1  = u + v;
        y2r = -(u+v)/2.;
        y2i =  (u-v)*sqrt(3.)/2.;
    }
    
    
    temp1 = b/a/3.;
    y1 = y1-temp1;
    if (DD < 0.)
    {    y2 = y2-temp1;
         y3 = y3-temp1;
    }
    else
        y2r=y2r-temp1;
    
    
    if (DD < 0.)
    {
        roots[0] = y1;
        roots[1] = y2;
        roots[2] = y3;
        
    }
    else if (DD == 0.){
        roots[0] =  y1;
        roots[1] = y2r;
        roots[2] = y2r;
    }
    else
    {
        roots[0] = y1;
        roots[1] = y2r;
        roots[2] = y2r;
        if (y2i != 0)
        {
            roots[1] = 0;
            roots[2] = 0;
        }
        
    }
    
}


/* main function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    
    srand((unsigned)time(NULL));
    int n;                          // size of the input
    double *D;                      // input matrix
    double *D_estim;                // estimated EDM
    double *dist;                   // distance to the estimated EDM
    double *XY_coor_estim_output;   // the estimated positions for the output
    int iter_max = 150;             // maximum number of iterations
    double** XY_coor_init;          // initial point for starting the minimization
    double** XY_coor_estim;         // estimated configuration
    double** D_mat;                 // matrix form of D;
    int iter_count;                 // itration counter
    int sen_ind;                    // sensor index counter
    double a, b, c, d;              // third order polynomial coeffs
    double cbcroots[3];             // third order polynomial roots
    double deltaX_min;              // optimization step in X direction.
    double deltaY_min;              // optimization step in Y direction.
    double deltaZ_min;              // optimization step in Z direction.
    int dim = 2;                    // Default value for the dimension
    double cost[4];
    double min_cost;
    int i;
    int j;
    int k;
    
    cost[0] = 10000;                // dummy value for finding the min cost
    
    
    /* define the input and output arrays */
    if (nrhs < 1 || nrhs > 2) 
    {
        mexErrMsgTxt("At least one and at most two input arguments required.");
    } 
    else if (nlhs > 3) 
    {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    if (nrhs > 1) 
    {
        dim = *mxGetPr(prhs[1]);
    }
    
    if (!(dim==2 || dim==3))
    {
        mexErrMsgTxt("This program supports only dimensions 2 or 3.");
    }
    
    
    
    D = mxGetPr(prhs[0]);           // input matrix
    n = mxGetM(prhs[0]);            // input dimension
    
    /* create pointee to dim */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    dist = mxGetPr(plhs[0]);
    
    /* create pointer to estimated EDM */
    plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
    D_estim = mxGetPr(plhs[1]);
    
    /* create pointer to estimated coordinates */
    plhs[2] = mxCreateDoubleMatrix(n, dim, mxREAL); // estimated XY coordinates
    XY_coor_estim_output = mxGetPr(plhs[2]);
    
    
    /* initialize the pointers for initial solution, estimated coordinates
     * and matrix form of the input matrix */
    XY_coor_init = (double**)mxMalloc(n * sizeof(double*));
    for (i = 0; i < n; i++) 
    {
        XY_coor_init[i] = (double*)mxMalloc(dim * sizeof(double));
    }
    
    XY_coor_estim = (double**)mxMalloc(n * sizeof(double*));
    for (i = 0; i < n; i++) 
    {
        XY_coor_estim[i] = (double*)mxMalloc(dim * sizeof(double));
    }
    
    D_mat = (double**)mxMalloc(n * sizeof(double*));
    for (i = 0; i < n; i++) 
    {
        D_mat[i] = (double*)mxMalloc(n * sizeof(double));
    }
    
    
    
    /* Initialize the starting point of the algorithm */
    for (i = 0; i < n; i++)
        for (j = 0; j < dim; j++)
        {
            XY_coor_init[i][j] = 0;
            XY_coor_estim[i][j] = XY_coor_init[i][j] ;
        }
    
    /* Transform the input to the matrix format */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            D_mat[i][j] = D[i*n+j];
    
    /* If dim===2, perform a 2-dimensional optimzation */
    if (2 == dim)
    {
        for (iter_count = 0; iter_count < iter_max; iter_count++)
        {
            for(sen_ind = 0; sen_ind < n; sen_ind++)
            {
                /* Update x coordinate */
                a = 4 * n;
                
                b = 0;
                c = 0;
                d = 0;

                for (j = 0; j < n; j++)
                {
                    b = b + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]);
                    c = c + 3*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) - D_mat[sen_ind][j];                 
                    d = d + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) * ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) - D_mat[sen_ind][j]);                
                }
                
                b = 12 * b;
                c = 4 * c;
                d = 4 * d;
                
                cubicroots(a, b, c, d, cbcroots);
                deltaX_min = cbcroots[0];
                min_cost = cost[0];
                for (k = 1; k < 4; k++)
                {
                    cost[k] = 0;
                    for (j = 0; j < n; j++)
                        if (j != sen_ind)
                            cost[k] = cost[k] + ((XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) - D_mat[sen_ind][j]) * ((XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) - D_mat[sen_ind][j]);
                    if (cost[k] < min_cost)
                    {
                        deltaX_min = cbcroots[k-1];
                        min_cost = cost[k];
                    }
                }
                XY_coor_estim[sen_ind][0] = XY_coor_estim[sen_ind][0] + deltaX_min;
                
                
                /* Update y coordinate */
                a = 4 * n;
                
                b = 0;
                for (j = 0; j < n; j++)
                    b = b + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]);
                b = 12 * b;
                
                c = 0;
                for (j = 0; j < n; j++)
                    c = c + 3*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) - D_mat[sen_ind][j];
                c = 4 * c;
                
                d = 0;
                for (j = 0; j < n; j++)
                    d = d + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) * ((XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])+(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])-D_mat[sen_ind][j]);
                d = 4 * d;
                
                cubicroots(a, b, c, d, cbcroots);
                deltaY_min = cbcroots[0];
                min_cost = cost[0];
                for (k = 1; k < 4; k++)
                {
                    cost[k] = 0;
                        for (j = 0; j < n; j++)
                            if (j != sen_ind)
                                cost[k] = cost[k] + ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1]+cbcroots[k-1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] +cbcroots[k-1] - XY_coor_estim[j][1]) - D_mat[sen_ind][j]) * ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1]+cbcroots[k-1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] +cbcroots[k-1] - XY_coor_estim[j][1]) - D_mat[sen_ind][j]);
                    if (cost[k] < min_cost)
                    {
                        deltaY_min = cbcroots[k-1];
                        min_cost = cost[k];
                    }
                }
                XY_coor_estim[sen_ind][1] = XY_coor_estim[sen_ind][1] + deltaY_min;
            } // endfor(sen_ind)
        } // endfor(iter_count)
    }
    /* If dim == 3 then perform a 3-dimensional optimization */
    else
    {
        for (iter_count = 0; iter_count < iter_max; iter_count++)
        {
            for(sen_ind = 0; sen_ind < n; sen_ind++)
            {
                
                /* Update x coordinate */
                a = 4 * n;
                
                b = 0;
                c = 0;
                d = 0;
                for (j = 0; j < n; j++)
                {
                    b = b + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]);
                    c = c + 3*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) * (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) - D_mat[sen_ind][j];
                    d = d + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) * ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) - D_mat[sen_ind][j]);
                }

                b = 12 * b;
                c = 4 * c;
                d = 4 * d;
                
                cubicroots(a, b, c, d, cbcroots);
                
                deltaX_min = cbcroots[0];
                min_cost = cost[0];
                for (k = 1; k < 4; k++)
                {
                    cost[k] = 0;
                    for (j = 0; j < n; j++)
                        if (j != sen_ind)
                            cost[k] = cost[k] + ((XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) - D_mat[sen_ind][j]) * ((XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0]) * (XY_coor_estim[sen_ind][0] + cbcroots[k-1] - XY_coor_estim[j][0]) + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) * (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) - D_mat[sen_ind][j]);
                    if (cost[k] < min_cost)
                    {
                        deltaX_min = cbcroots[k-1];
                        min_cost = cost[k];
                    }
                }
                XY_coor_estim[sen_ind][0] = XY_coor_estim[sen_ind][0] + deltaX_min;
                
                
                /* Update y coordinate */
                a = 4 * n;
                
                b = c = d = 0;
                
                for (j = 0; j < n; j++)
                {
                    b = b + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]);
                    c = c + 3*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])- D_mat[sen_ind][j];
                    d = d + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) * ((XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])+(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])-D_mat[sen_ind][j]);
                }

                b = 12 * b;
                c = 4 * c;
                d = 4 * d;
                
                cubicroots(a, b, c, d, cbcroots);
                deltaY_min = cbcroots[0];
                
                min_cost = cost[0];
                for (k = 1; k < 4; k++)
                {
                    cost[k] = 0;
                        for (j = 0; j < n; j++)
                            if (j != sen_ind)
                                cost[k] = cost[k] + ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1]+cbcroots[k-1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] +cbcroots[k-1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) - D_mat[sen_ind][j]) * ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1]+cbcroots[k-1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] +cbcroots[k-1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) - D_mat[sen_ind][j]);
                    if (cost[k] < min_cost)
                    {
                        deltaY_min = cbcroots[k-1];
                        min_cost = cost[k];
                    }
                }
                XY_coor_estim[sen_ind][1] = XY_coor_estim[sen_ind][1] + deltaY_min;
                
                /* Update z coordinate */
                a = 4 * n;
                
                b = c = d = 0;
                for (j = 0; j < n; j++)
                {
                    b = b + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]);
                    c = c + 3*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) + (XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])- D_mat[sen_ind][j];
                    d = d + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2]) * ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2] - XY_coor_estim[j][2])- D_mat[sen_ind][j]);
                }

                b = 12 * b;
                c = 4 * c;
                d = 4 * d;
                
                cubicroots(a, b, c, d, cbcroots);
                
                deltaZ_min = cbcroots[0];
                min_cost = cost[0];
                for (k = 1; k < 4; k++)
                {
                    cost[k] = 0;
                    for (j = 0; j < n; j++)
                        if (j != sen_ind)
                            cost[k] = cost[k] + ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2]  + cbcroots[k-1]- XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2]  + cbcroots[k-1]- XY_coor_estim[j][2])- D_mat[sen_ind][j]) * ((XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])*(XY_coor_estim[sen_ind][0] - XY_coor_estim[j][0])+(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1])*(XY_coor_estim[sen_ind][1] - XY_coor_estim[j][1]) + (XY_coor_estim[sen_ind][2]  + cbcroots[k-1]- XY_coor_estim[j][2])*(XY_coor_estim[sen_ind][2]  + cbcroots[k-1]- XY_coor_estim[j][2])- D_mat[sen_ind][j]);
                    if (cost[k] < min_cost)
                    {
                        deltaZ_min = cbcroots[k-1];
                        min_cost = cost[k];
                    }
                }
                XY_coor_estim[sen_ind][2] = XY_coor_estim[sen_ind][2] + deltaZ_min;
            } // endfor(sen_ind)
        } // endfor(iter_count)
    } // endelse
    
    
    /* Find the optimal distance matrix and the Frobenius of the error */
    if (dim == 2)
    {
        dist[0] = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                dist[0] = dist[0] + ((XY_coor_estim[i][0] - XY_coor_estim[j][0])*(XY_coor_estim[i][0] - XY_coor_estim[j][0]) + (XY_coor_estim[i][1] - XY_coor_estim[j][1])*(XY_coor_estim[i][1] - XY_coor_estim[j][1]) - D_mat[i][j])*((XY_coor_estim[i][0] - XY_coor_estim[j][0])*(XY_coor_estim[i][0] - XY_coor_estim[j][0]) + (XY_coor_estim[i][1] - XY_coor_estim[j][1])*(XY_coor_estim[i][1] - XY_coor_estim[j][1]) - D_mat[i][j]);
        
        dist[0] = 1./n * sqrt(dist[0]);
        
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                D_estim[i*n+j] = (XY_coor_estim[i][0] - XY_coor_estim[j][0])*(XY_coor_estim[i][0] - XY_coor_estim[j][0]) + (XY_coor_estim[i][1] - XY_coor_estim[j][1])*(XY_coor_estim[i][1] - XY_coor_estim[j][1]);
    }
    else
    {
        dist[0] = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                dist[0] = dist[0] + ((XY_coor_estim[i][0] - XY_coor_estim[j][0])*(XY_coor_estim[i][0] - XY_coor_estim[j][0]) + (XY_coor_estim[i][1] - XY_coor_estim[j][1])*(XY_coor_estim[i][1] - XY_coor_estim[j][1]) + (XY_coor_estim[i][2] - XY_coor_estim[j][2])*(XY_coor_estim[i][2] - XY_coor_estim[j][2])- D_mat[i][j])*((XY_coor_estim[i][0] - XY_coor_estim[j][0])*(XY_coor_estim[i][0] - XY_coor_estim[j][0]) + (XY_coor_estim[i][1] - XY_coor_estim[j][1])*(XY_coor_estim[i][1] - XY_coor_estim[j][1]) + (XY_coor_estim[i][2] - XY_coor_estim[j][2])*(XY_coor_estim[i][2] - XY_coor_estim[j][2])- D_mat[i][j]);
        
        dist[0] = 1./n * sqrt(dist[0]);
        
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                D_estim[i*n+j] = (XY_coor_estim[i][0] - XY_coor_estim[j][0])*(XY_coor_estim[i][0] - XY_coor_estim[j][0]) + (XY_coor_estim[i][1] - XY_coor_estim[j][1])*(XY_coor_estim[i][1] - XY_coor_estim[j][1]) + (XY_coor_estim[i][2] - XY_coor_estim[j][2])*(XY_coor_estim[i][2] - XY_coor_estim[j][2]);
        
    }
    for (j = 0; j < dim; j++)
        for (i = 0; i < n; i++)
            XY_coor_estim_output[i+j*n] = XY_coor_estim[i][j];
    
} // endfunction