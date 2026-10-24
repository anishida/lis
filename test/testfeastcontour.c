#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <stdio.h>
#include <math.h>

#include "lis.h"
#include "lis_esolver.h"


static int close_real(
        LIS_REAL a,
        LIS_REAL b,
        LIS_REAL tol)
{
        LIS_REAL scale;

        scale = (LIS_REAL)1.0+fabs(b);

        return fabs(a-b) <= tol*scale;
}


static int test_gauss_legendre(void)
{
        LIS_REAL x[8];
        LIS_REAL w[8];

        LIS_REAL m0;
        LIS_REAL m1;
        LIS_REAL m2;
        LIS_REAL m14;

        LIS_INT i;
        LIS_INT err;

        m0 = 0.0;
        m1 = 0.0;
        m2 = 0.0;
        m14 = 0.0;

        err = lis_efeast_gauss_legendre(
                8,
                x,
                w);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: Gauss-Legendre generation returned %d\n",
                        (int)err);

                return 1;
        }

        for(i=0;i<8;i++)
        {
                LIS_REAL x2;
                LIS_REAL x4;
                LIS_REAL x8;
                LIS_REAL x14;

                if( w[i]<=0.0 )
                {
                        fprintf(
                                stderr,
                                "FAIL: non-positive Gauss weight at %d\n",
                                (int)i);

                        return 1;
                }

                if( i>0 && !(x[i-1]<x[i]) )
                {
                        fprintf(
                                stderr,
                                "FAIL: Gauss nodes are not strictly ordered\n");

                        return 1;
                }

                x2 = x[i]*x[i];
                x4 = x2*x2;
                x8 = x4*x4;
                x14 = x8*x4*x2;

                m0 += w[i];
                m1 += w[i]*x[i];
                m2 += w[i]*x2;
                m14 += w[i]*x14;
        }

        if( !close_real(
                    m0,
                    (LIS_REAL)2.0,
                    (LIS_REAL)2.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral 1 = %.17g\n",
                        (double)m0);

                return 1;
        }

        if( !close_real(
                    m1,
                    (LIS_REAL)0.0,
                    (LIS_REAL)2.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral x = %.17g\n",
                        (double)m1);

                return 1;
        }

        if( !close_real(
                    m2,
                    (LIS_REAL)(2.0/3.0),
                    (LIS_REAL)2.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral x^2 = %.17g\n",
                        (double)m2);

                return 1;
        }

        if( !close_real(
                    m14,
                    (LIS_REAL)(2.0/15.0),
                    (LIS_REAL)3.0e-13) )
        {
                fprintf(
                        stderr,
                        "FAIL: integral x^14 = %.17g\n",
                        (double)m14);

                return 1;
        }

        return 0;
}


static int test_feast_reference_contour(void)
{
        LIS_REAL zr[4];
        LIS_REAL zi[4];
        LIS_REAL wr[4];
        LIS_REAL wi[4];

        const LIS_REAL ref_zr[4] = {
                (LIS_REAL)1.0473908518410222,
                (LIS_REAL)1.9819684290084905,
                (LIS_REAL)4.0180315709915098,
                (LIS_REAL)4.9526091481589782
        };

        const LIS_REAL ref_zi[4] = {
                (LIS_REAL)0.12984058035656038,
                (LIS_REAL)0.51645431050753199,
                (LIS_REAL)0.51645431050753199,
                (LIS_REAL)0.12984058035656015
        };

        const LIS_REAL ref_wr[4] = {
                (LIS_REAL)-0.050941841463511224,
                (LIS_REAL)-0.049792826738941448,
                (LIS_REAL) 0.049792826738941448,
                (LIS_REAL) 0.050941841463511231
        };

        const LIS_REAL ref_wi[4] = {
                (LIS_REAL)0.037638062477073672,
                (LIS_REAL)0.28066931358780345,
                (LIS_REAL)0.28066931358780339,
                (LIS_REAL)0.037638062477073610
        };

        LIS_INT i;
        LIS_INT err;

        err = lis_efeast_contour_gauss(
                (LIS_REAL)1.0,
                (LIS_REAL)5.0,
                4,
                (LIS_REAL)0.30,
                zr,
                zi,
                wr,
                wi);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: FEAST contour generation returned %d\n",
                        (int)err);

                return 1;
        }

        for(i=0;i<4;i++)
        {
                if( !close_real(
                            zr[i],
                            ref_zr[i],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            zi[i],
                            ref_zi[i],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wr[i],
                            ref_wr[i],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wi[i],
                            ref_wi[i],
                            (LIS_REAL)5.0e-13) )
                {
                        fprintf(
                                stderr,
                                "FAIL: contour node %d differs from reference\n"
                                "  z=(%.17g, %.17g)\n"
                                "  w=(%.17g, %.17g)\n",
                                (int)i,
                                (double)zr[i],
                                (double)zi[i],
                                (double)wr[i],
                                (double)wi[i]);

                        return 1;
                }
        }

        return 0;
}


static int test_contour_symmetry(void)
{
        LIS_REAL zr[8];
        LIS_REAL zi[8];
        LIS_REAL wr[8];
        LIS_REAL wi[8];

        LIS_INT i;
        LIS_INT j;
        LIS_INT err;

        err = lis_efeast_contour_gauss(
                (LIS_REAL)-3.0,
                (LIS_REAL)7.0,
                8,
                (LIS_REAL)0.30,
                zr,
                zi,
                wr,
                wi);

        if( err )
        {
                return 1;
        }

        for(i=0;i<4;i++)
        {
                j = 7-i;

                if( !close_real(
                            zr[i]+zr[j],
                            (LIS_REAL)4.0,
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            zi[i],
                            zi[j],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wr[i],
                            -wr[j],
                            (LIS_REAL)5.0e-13) ||
                    !close_real(
                            wi[i],
                            wi[j],
                            (LIS_REAL)5.0e-13) )
                {
                        fprintf(
                                stderr,
                                "FAIL: contour symmetry at pair %d/%d\n",
                                (int)i,
                                (int)j);

                        return 1;
                }
        }

        return 0;
}


#define PROJECTOR_TEST_N 4


static LIS_REAL projector_scalar_real(
        LIS_SCALAR value)
{
#ifdef _COMPLEX
        return (LIS_REAL)creal(value);
#else
        return (LIS_REAL)value;
#endif
}


static LIS_REAL projector_scalar_imag(
        LIS_SCALAR value)
{
#ifdef _COMPLEX
        return (LIS_REAL)cimag(value);
#else
        (void)value;
        return (LIS_REAL)0.0;
#endif
}


static LIS_INT create_projector_diagonal(
        const LIS_REAL diagonal[PROJECTOR_TEST_N],
        LIS_MATRIX *matrix)
{
        LIS_MATRIX A;

        LIS_INT i;
        LIS_INT err;

        A = NULL;

        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &A);

        if( err ) goto fail;

        err = lis_matrix_set_size(
                A,
                0,
                PROJECTOR_TEST_N);

        if( err ) goto fail;

        for(i=0;i<PROJECTOR_TEST_N;i++)
        {
                err = lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        i,
                        (LIS_SCALAR)diagonal[i],
                        A);

                if( err ) goto fail;
        }

        err = lis_matrix_assemble(
                A);

        if( err ) goto fail;

        *matrix = A;

        return LIS_SUCCESS;

fail:

        if( A!=NULL )
        {
                lis_matrix_destroy(
                        A);
        }

        return err;
}


static LIS_INT set_projector_vector(
        LIS_VECTOR vector,
        const LIS_REAL values[PROJECTOR_TEST_N])
{
        LIS_INT i;
        LIS_INT err;

        for(i=0;i<PROJECTOR_TEST_N;i++)
        {
                err = lis_vector_set_value(
                        LIS_INS_VALUE,
                        i,
                        (LIS_SCALAR)values[i],
                        vector);

                if( err )
                {
                        return err;
                }
        }

        return LIS_SUCCESS;
}


static int check_projector_vector(
        const char *name,
        LIS_VECTOR vector,
        const LIS_REAL expected[PROJECTOR_TEST_N],
        LIS_REAL tol)
{
        LIS_SCALAR value;

        LIS_REAL actual_real;
        LIS_REAL actual_imag;

        LIS_INT i;
        LIS_INT err;

        for(i=0;i<PROJECTOR_TEST_N;i++)
        {
                err = lis_vector_get_value(
                        vector,
                        i,
                        &value);

                if( err )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s get value %d returned %d\n",
                                name,
                                (int)i,
                                (int)err);

                        return 1;
                }

                actual_real =
                        projector_scalar_real(
                                value);

                actual_imag =
                        projector_scalar_imag(
                                value);

                if( !close_real(
                            actual_real,
                            expected[i],
                            tol) ||
                    fabs(actual_imag) >
                            tol*((LIS_REAL)1.0+
                                 fabs(expected[i])) )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s[%d] = "
                                "(%.17g, %.17g), expected %.17g\n",
                                name,
                                (int)i,
                                (double)actual_real,
                                (double)actual_imag,
                                (double)expected[i]);

                        return 1;
                }
        }

        return 0;
}


static int run_projector_case(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B)
{
        const LIS_REAL input_real[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)3.0,
                (LIS_REAL)4.0
        };

        const LIS_REAL input_imag[PROJECTOR_TEST_N] = {
                (LIS_REAL)-0.5,
                (LIS_REAL)1.0,
                (LIS_REAL)-2.0,
                (LIS_REAL)0.25
        };

        const LIS_REAL expected_real[PROJECTOR_TEST_N] = {
                (LIS_REAL)0.0,
                (LIS_REAL)2.0,
                (LIS_REAL)3.0,
                (LIS_REAL)0.0
        };

        const LIS_REAL expected_imag[PROJECTOR_TEST_N] = {
                (LIS_REAL)0.0,
                (LIS_REAL)1.0,
                (LIS_REAL)-2.0,
                (LIS_REAL)0.0
        };

        LIS_VECTOR x_real;
        LIS_VECTOR x_imag;

        LIS_VECTOR q_real;
        LIS_VECTOR q_imag;

        LIS_SOLVER solver;

        LIS_INT err;

        x_real = NULL;
        x_imag = NULL;

        q_real = NULL;
        q_imag = NULL;

        solver = NULL;

        err = lis_vector_duplicate(
                A,
                &x_real);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &x_imag);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &q_real);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &q_imag);

        if( err ) goto fail;

        err = set_projector_vector(
                x_real,
                input_real);

        if( err ) goto fail;

        err = set_projector_vector(
                x_imag,
                input_imag);

        if( err ) goto fail;

        err = lis_solver_create(
                &solver);

        if( err ) goto fail;

        solver->options[
                LIS_OPTIONS_SOLVER] =
                LIS_SOLVER_GMRES;

        solver->options[
                LIS_OPTIONS_PRECON] =
                LIS_PRECON_TYPE_NONE;

        solver->options[
                LIS_OPTIONS_SCALE] =
                LIS_SCALE_NONE;

        solver->options[
                LIS_OPTIONS_STORAGE] =
                0;

        solver->options[
                LIS_OPTIONS_USE_AT] =
                LIS_FALSE;

        solver->options[
                LIS_OPTIONS_PRECISION] =
                LIS_PRECISION_DOUBLE;

        solver->options[
                LIS_OPTIONS_MAXITER] =
                200;

        solver->options[
                LIS_OPTIONS_RESTART] =
                16;

        solver->options[
                LIS_OPTIONS_INITGUESS_ZEROS] =
                LIS_TRUE;

        solver->params[
                LIS_PARAMS_RESID
                -
                LIS_OPTIONS_LEN] =
                (LIS_SCALAR)1.0e-13;

        /*
         * Eigenvalues are 1, 2, 4 and 7.
         *
         * The contour [1.5,4.5] therefore selects exactly
         * eigenvalues 2 and 4.
         *
         * Thirty-two upper-half Gauss points make the projector
         * error comfortably smaller than the regression tolerance.
         */
        err = lis_efeast_projector_apply_split(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)0.30,
                32,
                x_real,
                x_imag,
                q_real,
                q_imag,
                solver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s projector returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        if( check_projector_vector(
                    "projector real",
                    q_real,
                    expected_real,
                    (LIS_REAL)5.0e-6) )
        {
                goto fail;
        }

        if( check_projector_vector(
                    "projector imag",
                    q_imag,
                    expected_imag,
                    (LIS_REAL)5.0e-6) )
        {
                goto fail;
        }

        printf(
                "PASS: %s spectral projector\n",
                name);

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        if( q_imag!=NULL )
        {
                lis_vector_destroy(
                        q_imag);
        }

        if( q_real!=NULL )
        {
                lis_vector_destroy(
                        q_real);
        }

        if( x_imag!=NULL )
        {
                lis_vector_destroy(
                        x_imag);
        }

        if( x_real!=NULL )
        {
                lis_vector_destroy(
                        x_real);
        }

        return 0;

fail:

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        if( q_imag!=NULL )
        {
                lis_vector_destroy(
                        q_imag);
        }

        if( q_real!=NULL )
        {
                lis_vector_destroy(
                        q_real);
        }

        if( x_imag!=NULL )
        {
                lis_vector_destroy(
                        x_imag);
        }

        if( x_real!=NULL )
        {
                lis_vector_destroy(
                        x_real);
        }

        return 1;
}


static int test_spectral_projector(void)
{
        const LIS_REAL standard_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)4.0,
                (LIS_REAL)7.0
        };

        /*
         * Generalized A x = lambda B x with the same eigenvalues:
         *
         * A = diag(1,4,12,28)
         * B = diag(1,2,3,4)
         */
        const LIS_REAL generalized_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)4.0,
                (LIS_REAL)12.0,
                (LIS_REAL)28.0
        };

        const LIS_REAL generalized_b[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)3.0,
                (LIS_REAL)4.0
        };

        LIS_MATRIX A_standard;
        LIS_MATRIX A_generalized;
        LIS_MATRIX B_generalized;

        LIS_INT err;

        A_standard = NULL;
        A_generalized = NULL;
        B_generalized = NULL;

        err = create_projector_diagonal(
                standard_a,
                &A_standard);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_a,
                &A_generalized);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_b,
                &B_generalized);

        if( err ) goto fail;

        if( run_projector_case(
                    "standard",
                    A_standard,
                    NULL) )
        {
                goto fail;
        }

        if( run_projector_case(
                    "generalized",
                    A_generalized,
                    B_generalized) )
        {
                goto fail;
        }

        lis_matrix_destroy(
                B_generalized);

        lis_matrix_destroy(
                A_generalized);

        lis_matrix_destroy(
                A_standard);

        printf(
                "PASS: native FEAST spectral projector "
                "standard/generalized\n");

        return 0;

fail:

        if( B_generalized!=NULL )
        {
                lis_matrix_destroy(
                        B_generalized);
        }

        if( A_generalized!=NULL )
        {
                lis_matrix_destroy(
                        A_generalized);
        }

        if( A_standard!=NULL )
        {
                lis_matrix_destroy(
                        A_standard);
        }

        return 1;
}


#define FILTERED_BLOCK_SIZE 3


static LIS_INT set_filtered_candidate(
        LIS_VECTOR vector,
        LIS_INT candidate)
{
        LIS_INT err;

        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                vector);

        if( err )
        {
                return err;
        }

        switch(candidate)
        {
        case 0:

#ifdef _COMPLEX

                return lis_vector_set_value(
                        LIS_INS_VALUE,
                        1,
                        (LIS_SCALAR)(1.0 + 1.0*I),
                        vector);

#else

                return lis_vector_set_value(
                        LIS_INS_VALUE,
                        1,
                        (LIS_SCALAR)1.0,
                        vector);

#endif

        case 1:

#ifdef _COMPLEX

                return lis_vector_set_value(
                        LIS_INS_VALUE,
                        2,
                        (LIS_SCALAR)(2.0 - 1.0*I),
                        vector);

#else

                return lis_vector_set_value(
                        LIS_INS_VALUE,
                        2,
                        (LIS_SCALAR)2.0,
                        vector);

#endif

        case 2:

#ifdef _COMPLEX

                err = lis_vector_set_value(
                        LIS_INS_VALUE,
                        1,
                        (LIS_SCALAR)(1.0 + 1.0*I),
                        vector);

                if( err )
                {
                        return err;
                }

                return lis_vector_set_value(
                        LIS_INS_VALUE,
                        2,
                        (LIS_SCALAR)(2.0 - 1.0*I),
                        vector);

#else

                err = lis_vector_set_value(
                        LIS_INS_VALUE,
                        1,
                        (LIS_SCALAR)1.0,
                        vector);

                if( err )
                {
                        return err;
                }

                return lis_vector_set_value(
                        LIS_INS_VALUE,
                        2,
                        (LIS_SCALAR)2.0,
                        vector);

#endif

        default:

                return LIS_ERR_ILL_ARG;
        }
}


static int check_filtered_metric_identity(
        LIS_INT rank,
        LIS_VECTOR Q[],
        LIS_VECTOR BQ[],
        LIS_REAL tol)
{
        LIS_SCALAR value;

        LIS_REAL expected;
        LIS_REAL real_part;
        LIS_REAL imag_part;

        LIS_INT i;
        LIS_INT j;
        LIS_INT err;

        for(i=0;i<rank;i++)
        {
                for(j=0;j<rank;j++)
                {
                        if( BQ!=NULL )
                        {
                                err = lis_vector_dot(
                                        Q[i],
                                        BQ[j],
                                        &value);
                        }
                        else
                        {
                                err = lis_vector_dot(
                                        Q[i],
                                        Q[j],
                                        &value);
                        }

                        if( err )
                        {
                                fprintf(
                                        stderr,
                                        "FAIL: filtered metric dot returned %d\n",
                                        (int)err);

                                return 1;
                        }

                        expected =
                                (i==j)
                                ?
                                (LIS_REAL)1.0
                                :
                                (LIS_REAL)0.0;

                        real_part =
                                projector_scalar_real(
                                        value);

                        imag_part =
                                projector_scalar_imag(
                                        value);

                        if( fabs(real_part-expected)>tol ||
                            fabs(imag_part)>tol )
                        {
                                fprintf(
                                        stderr,
                                        "FAIL: filtered metric (%d,%d) = "
                                        "(%.17g, %.17g), expected %.17g\n",
                                        (int)i,
                                        (int)j,
                                        (double)real_part,
                                        (double)imag_part,
                                        (double)expected);

                                return 1;
                        }
                }
        }

        return 0;
}


static int check_filtered_support(
        LIS_INT rank,
        LIS_VECTOR Q[],
        LIS_REAL tol)
{
        LIS_SCALAR value;

        LIS_REAL magnitude;

        LIS_INT i;
        LIS_INT index;
        LIS_INT err;

        const LIS_INT outside[2] = {
                0,
                3
        };

        for(i=0;i<rank;i++)
        {
                for(index=0;index<2;index++)
                {
                        err = lis_vector_get_value(
                                Q[i],
                                outside[index],
                                &value);

                        if( err )
                        {
                                fprintf(
                                        stderr,
                                        "FAIL: filtered support read returned %d\n",
                                        (int)err);

                                return 1;
                        }

                        magnitude = sqrt(
                                projector_scalar_real(value)
                                *
                                projector_scalar_real(value)
                                +
                                projector_scalar_imag(value)
                                *
                                projector_scalar_imag(value));

                        if( magnitude>tol )
                        {
                                fprintf(
                                        stderr,
                                        "FAIL: filtered basis leaked to "
                                        "component %d: %.17g\n",
                                        (int)outside[index],
                                        (double)magnitude);

                                return 1;
                        }
                }
        }

        return 0;
}


static int check_filtered_zero_tail(
        LIS_VECTOR vector,
        LIS_REAL tol)
{
        LIS_REAL norm;

        LIS_INT err;

        norm = 0.0;

        err = lis_vector_nrm2(
                vector,
                &norm);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: filtered zero-tail norm returned %d\n",
                        (int)err);

                return 1;
        }

        if( norm>tol )
        {
                fprintf(
                        stderr,
                        "FAIL: filtered inactive vector norm %.17g\n",
                        (double)norm);

                return 1;
        }

        return 0;
}


static int run_filtered_block_case(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B)
{
        LIS_VECTOR X[FILTERED_BLOCK_SIZE];
        LIS_VECTOR Q[FILTERED_BLOCK_SIZE];
        LIS_VECTOR BQ[FILTERED_BLOCK_SIZE];

        LIS_SOLVER solver;

        LIS_INT active_rank;
        LIS_INT i;
        LIS_INT err;

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                X[i] = NULL;
                Q[i] = NULL;
                BQ[i] = NULL;
        }

        solver = NULL;
        active_rank = -1;

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &X[i]);

                if( err ) goto fail;

                err = lis_vector_duplicate(
                        A,
                        &Q[i]);

                if( err ) goto fail;

                if( B!=NULL )
                {
                        err = lis_vector_duplicate(
                                A,
                                &BQ[i]);

                        if( err ) goto fail;
                }

                err = set_filtered_candidate(
                        X[i],
                        i);

                if( err ) goto fail;
        }

        err = lis_solver_create(
                &solver);

        if( err ) goto fail;

        solver->options[
                LIS_OPTIONS_SOLVER] =
                LIS_SOLVER_GMRES;

        solver->options[
                LIS_OPTIONS_PRECON] =
                LIS_PRECON_TYPE_NONE;

        solver->options[
                LIS_OPTIONS_SCALE] =
                LIS_SCALE_NONE;

        solver->options[
                LIS_OPTIONS_STORAGE] =
                0;

        solver->options[
                LIS_OPTIONS_USE_AT] =
                LIS_FALSE;

        solver->options[
                LIS_OPTIONS_PRECISION] =
                LIS_PRECISION_DOUBLE;

        solver->options[
                LIS_OPTIONS_MAXITER] =
                200;

        solver->options[
                LIS_OPTIONS_RESTART] =
                16;

        solver->options[
                LIS_OPTIONS_INITGUESS_ZEROS] =
                LIS_TRUE;

        solver->params[
                LIS_PARAMS_RESID
                -
                LIS_OPTIONS_LEN] =
                (LIS_SCALAR)1.0e-13;

        /*
         * The projected candidates span exactly the two-dimensional
         * invariant subspace associated with eigenvalues 2 and 4:
         *
         *     X0 -> e_1
         *     X1 -> e_2
         *     X2 -> X0 + X1
         *
         * Therefore active_rank must be two.
         */
        err = lis_efeast_filtered_block_basis(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)0.30,
                32,
                FILTERED_BLOCK_SIZE,
                X,
                Q,
                B!=NULL ? BQ : NULL,
                (LIS_REAL)1.0e-8,
                &active_rank,
                solver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s filtered block returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        if( active_rank!=2 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s active_rank=%d expected 2\n",
                        name,
                        (int)active_rank);

                goto fail;
        }

        if( check_filtered_metric_identity(
                    active_rank,
                    Q,
                    B!=NULL ? BQ : NULL,
                    (LIS_REAL)2.0e-8) )
        {
                goto fail;
        }

        if( check_filtered_support(
                    active_rank,
                    Q,
                    (LIS_REAL)2.0e-7) )
        {
                goto fail;
        }

        if( check_filtered_zero_tail(
                    Q[2],
                    (LIS_REAL)2.0e-10) )
        {
                goto fail;
        }

        if( B!=NULL &&
            check_filtered_zero_tail(
                    BQ[2],
                    (LIS_REAL)2.0e-10) )
        {
                goto fail;
        }

        printf(
                "PASS: %s filtered block basis active_rank=%d\n",
                name,
                (int)active_rank);

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                if( BQ[i]!=NULL )
                {
                        lis_vector_destroy(
                                BQ[i]);
                }

                if( Q[i]!=NULL )
                {
                        lis_vector_destroy(
                                Q[i]);
                }

                if( X[i]!=NULL )
                {
                        lis_vector_destroy(
                                X[i]);
                }
        }

        return 0;

fail:

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                if( BQ[i]!=NULL )
                {
                        lis_vector_destroy(
                                BQ[i]);
                }

                if( Q[i]!=NULL )
                {
                        lis_vector_destroy(
                                Q[i]);
                }

                if( X[i]!=NULL )
                {
                        lis_vector_destroy(
                                X[i]);
                }
        }

        return 1;
}


static int test_filtered_block_basis(void)
{
        const LIS_REAL standard_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)4.0,
                (LIS_REAL)7.0
        };

        const LIS_REAL generalized_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)4.0,
                (LIS_REAL)12.0,
                (LIS_REAL)28.0
        };

        const LIS_REAL generalized_b[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)3.0,
                (LIS_REAL)4.0
        };

        LIS_MATRIX A_standard;
        LIS_MATRIX A_generalized;
        LIS_MATRIX B_generalized;

        LIS_INT err;

        A_standard = NULL;
        A_generalized = NULL;
        B_generalized = NULL;

        err = create_projector_diagonal(
                standard_a,
                &A_standard);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_a,
                &A_generalized);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_b,
                &B_generalized);

        if( err ) goto fail;

        if( run_filtered_block_case(
                    "standard",
                    A_standard,
                    NULL) )
        {
                goto fail;
        }

        if( run_filtered_block_case(
                    "generalized",
                    A_generalized,
                    B_generalized) )
        {
                goto fail;
        }

        lis_matrix_destroy(
                B_generalized);

        lis_matrix_destroy(
                A_generalized);

        lis_matrix_destroy(
                A_standard);

        printf(
                "PASS: native FEAST filtered block basis "
                "standard/generalized active-rank\n");

        return 0;

fail:

        if( B_generalized!=NULL )
        {
                lis_matrix_destroy(
                        B_generalized);
        }

        if( A_generalized!=NULL )
        {
                lis_matrix_destroy(
                        A_generalized);
        }

        if( A_standard!=NULL )
        {
                lis_matrix_destroy(
                        A_standard);
        }

        return 1;
}


#define RR_TEST_RANK 2


static int check_rayleigh_ritz_values(
        const char *name,
        LIS_SCALAR values[RR_TEST_RANK],
        LIS_REAL tol)
{
        LIS_REAL a;
        LIS_REAL b;
        LIS_REAL ai;
        LIS_REAL bi;
        LIS_REAL low;
        LIS_REAL high;

        a = projector_scalar_real(
                values[0]);

        b = projector_scalar_real(
                values[1]);

        ai = projector_scalar_imag(
                values[0]);

        bi = projector_scalar_imag(
                values[1]);

        if( fabs(ai)>tol ||
            fabs(bi)>tol )
        {
                fprintf(
                        stderr,
                        "FAIL: %s projected eigenvalue "
                        "has non-real part\n",
                        name);

                return 1;
        }

        if( a<b )
        {
                low = a;
                high = b;
        }
        else
        {
                low = b;
                high = a;
        }

        if( fabs(low-(LIS_REAL)2.0)>tol ||
            fabs(high-(LIS_REAL)4.0)>tol )
        {
                fprintf(
                        stderr,
                        "FAIL: %s projected eigenvalues "
                        "%.17g %.17g, expected 2 and 4\n",
                        name,
                        (double)low,
                        (double)high);

                return 1;
        }

        return 0;
}


static int check_rayleigh_ritz_residuals(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_SCALAR values[RR_TEST_RANK],
        LIS_VECTOR vectors[RR_TEST_RANK],
        LIS_REAL tol)
{
        LIS_VECTOR Av;
        LIS_VECTOR Bv;

        LIS_REAL norm;

        LIS_INT i;
        LIS_INT err;

        Av = NULL;
        Bv = NULL;

        err = lis_vector_duplicate(
                A,
                &Av);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &Bv);

        if( err ) goto fail;

        for(i=0;i<RR_TEST_RANK;i++)
        {
                err = lis_matvec(
                        A,
                        vectors[i],
                        Av);

                if( err ) goto fail;

                if( B!=NULL )
                {
                        err = lis_matvec(
                                B,
                                vectors[i],
                                Bv);

                        if( err ) goto fail;
                }
                else
                {
                        err = lis_vector_copy(
                                vectors[i],
                                Bv);

                        if( err ) goto fail;
                }

                err = lis_vector_axpy(
                        -values[i],
                        Bv,
                        Av);

                if( err ) goto fail;

                norm = 0.0;

                err = lis_vector_nrm2(
                        Av,
                        &norm);

                if( err ) goto fail;

                if( !isfinite((double)norm) ||
                    norm>tol )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s Ritz residual %d = %.17g\n",
                                name,
                                (int)i,
                                (double)norm);

                        goto fail;
                }
        }

        lis_vector_destroy(
                Bv);

        lis_vector_destroy(
                Av);

        return 0;

fail:

        if( Bv!=NULL )
        {
                lis_vector_destroy(
                        Bv);
        }

        if( Av!=NULL )
        {
                lis_vector_destroy(
                        Av);
        }

        return 1;
}


static int check_rayleigh_ritz_metric(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR vectors[RR_TEST_RANK],
        LIS_REAL tol)
{
        LIS_VECTOR BV[RR_TEST_RANK];

        LIS_SCALAR dot;

        LIS_REAL expected;
        LIS_REAL real_part;
        LIS_REAL imag_part;

        LIS_INT i;
        LIS_INT j;
        LIS_INT err;

        for(i=0;i<RR_TEST_RANK;i++)
        {
                BV[i] = NULL;
        }

        if( B!=NULL )
        {
                for(i=0;i<RR_TEST_RANK;i++)
                {
                        err = lis_vector_duplicate(
                                A,
                                &BV[i]);

                        if( err ) goto fail;
                }

                err = lis_matvec_block(
                        B,
                        RR_TEST_RANK,
                        vectors,
                        BV);

                if( err ) goto fail;
        }

        for(i=0;i<RR_TEST_RANK;i++)
        {
                for(j=0;j<RR_TEST_RANK;j++)
                {
                        if( B!=NULL )
                        {
                                err = lis_vector_dot(
                                        vectors[i],
                                        BV[j],
                                        &dot);
                        }
                        else
                        {
                                err = lis_vector_dot(
                                        vectors[i],
                                        vectors[j],
                                        &dot);
                        }

                        if( err ) goto fail;

                        expected =
                                (i==j)
                                ?
                                (LIS_REAL)1.0
                                :
                                (LIS_REAL)0.0;

                        real_part =
                                projector_scalar_real(
                                        dot);

                        imag_part =
                                projector_scalar_imag(
                                        dot);

                        if( fabs(real_part-expected)>tol ||
                            fabs(imag_part)>tol )
                        {
                                fprintf(
                                        stderr,
                                        "FAIL: %s Ritz metric (%d,%d) = "
                                        "(%.17g, %.17g), expected %.17g\n",
                                        name,
                                        (int)i,
                                        (int)j,
                                        (double)real_part,
                                        (double)imag_part,
                                        (double)expected);

                                goto fail;
                        }
                }
        }

        for(i=0;i<RR_TEST_RANK;i++)
        {
                if( BV[i]!=NULL )
                {
                        lis_vector_destroy(
                                BV[i]);
                }
        }

        return 0;

fail:

        for(i=0;i<RR_TEST_RANK;i++)
        {
                if( BV[i]!=NULL )
                {
                        lis_vector_destroy(
                                BV[i]);
                }
        }

        return 1;
}


static int run_rayleigh_ritz_case(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B)
{
        LIS_VECTOR X[FILTERED_BLOCK_SIZE];
        LIS_VECTOR Q[FILTERED_BLOCK_SIZE];
        LIS_VECTOR BQ[FILTERED_BLOCK_SIZE];

        LIS_VECTOR ritz[RR_TEST_RANK];

        LIS_SCALAR values[RR_TEST_RANK];

        LIS_SOLVER solver;

        LIS_REAL projected_error;

        LIS_INT projected_niter;
        LIS_INT active_rank;
        LIS_INT i;
        LIS_INT err;

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                X[i] = NULL;
                Q[i] = NULL;
                BQ[i] = NULL;
        }

        for(i=0;i<RR_TEST_RANK;i++)
        {
                ritz[i] = NULL;
                values[i] = 0.0;
        }

        solver = NULL;

        active_rank = -1;
        projected_niter = 0;
        projected_error = 0.0;

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &X[i]);

                if( err ) goto fail;

                err = lis_vector_duplicate(
                        A,
                        &Q[i]);

                if( err ) goto fail;

                if( B!=NULL )
                {
                        err = lis_vector_duplicate(
                                A,
                                &BQ[i]);

                        if( err ) goto fail;
                }

                err = set_filtered_candidate(
                        X[i],
                        i);

                if( err ) goto fail;
        }

        for(i=0;i<RR_TEST_RANK;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &ritz[i]);

                if( err ) goto fail;
        }

        err = lis_solver_create(
                &solver);

        if( err ) goto fail;

        solver->options[
                LIS_OPTIONS_SOLVER] =
                LIS_SOLVER_GMRES;

        solver->options[
                LIS_OPTIONS_PRECON] =
                LIS_PRECON_TYPE_NONE;

        solver->options[
                LIS_OPTIONS_SCALE] =
                LIS_SCALE_NONE;

        solver->options[
                LIS_OPTIONS_STORAGE] =
                0;

        solver->options[
                LIS_OPTIONS_USE_AT] =
                LIS_FALSE;

        solver->options[
                LIS_OPTIONS_PRECISION] =
                LIS_PRECISION_DOUBLE;

        solver->options[
                LIS_OPTIONS_MAXITER] =
                200;

        solver->options[
                LIS_OPTIONS_RESTART] =
                16;

        solver->options[
                LIS_OPTIONS_INITGUESS_ZEROS] =
                LIS_TRUE;

        solver->params[
                LIS_PARAMS_RESID
                -
                LIS_OPTIONS_LEN] =
                (LIS_SCALAR)1.0e-13;

        /*
         * Build the same filtered two-dimensional invariant
         * subspace used by the F3-4 regression.
         */
        err = lis_efeast_filtered_block_basis(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)0.30,
                32,
                FILTERED_BLOCK_SIZE,
                X,
                Q,
                B!=NULL ? BQ : NULL,
                (LIS_REAL)1.0e-8,
                &active_rank,
                solver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s RR filtered basis returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        if( active_rank!=RR_TEST_RANK )
        {
                fprintf(
                        stderr,
                        "FAIL: %s RR active_rank=%d expected %d\n",
                        name,
                        (int)active_rank,
                        (int)RR_TEST_RANK);

                goto fail;
        }

        err = lis_efeast_projected_rayleigh_ritz(
                A,
                active_rank,
                Q,
                values,
                ritz,
                &projected_niter,
                &projected_error);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s projected Rayleigh-Ritz returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        if( projected_niter<=0 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s projected iteration count %d\n",
                        name,
                        (int)projected_niter);

                goto fail;
        }

        if( !isfinite((double)projected_error) )
        {
                fprintf(
                        stderr,
                        "FAIL: %s projected error is non-finite\n",
                        name);

                goto fail;
        }

        if( check_rayleigh_ritz_values(
                    name,
                    values,
                    (LIS_REAL)5.0e-6) )
        {
                goto fail;
        }

        if( check_rayleigh_ritz_residuals(
                    name,
                    A,
                    B,
                    values,
                    ritz,
                    (LIS_REAL)5.0e-6) )
        {
                goto fail;
        }

        if( check_rayleigh_ritz_metric(
                    name,
                    A,
                    B,
                    ritz,
                    (LIS_REAL)5.0e-7) )
        {
                goto fail;
        }

        printf(
                "PASS: %s projected Rayleigh-Ritz "
                "rank=%d iter=%d error=%.6e\n",
                name,
                (int)active_rank,
                (int)projected_niter,
                (double)projected_error);

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        for(i=0;i<RR_TEST_RANK;i++)
        {
                if( ritz[i]!=NULL )
                {
                        lis_vector_destroy(
                                ritz[i]);
                }
        }

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                if( BQ[i]!=NULL )
                {
                        lis_vector_destroy(
                                BQ[i]);
                }

                if( Q[i]!=NULL )
                {
                        lis_vector_destroy(
                                Q[i]);
                }

                if( X[i]!=NULL )
                {
                        lis_vector_destroy(
                                X[i]);
                }
        }

        return 0;

fail:

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        for(i=0;i<RR_TEST_RANK;i++)
        {
                if( ritz[i]!=NULL )
                {
                        lis_vector_destroy(
                                ritz[i]);
                }
        }

        for(i=0;i<FILTERED_BLOCK_SIZE;i++)
        {
                if( BQ[i]!=NULL )
                {
                        lis_vector_destroy(
                                BQ[i]);
                }

                if( Q[i]!=NULL )
                {
                        lis_vector_destroy(
                                Q[i]);
                }

                if( X[i]!=NULL )
                {
                        lis_vector_destroy(
                                X[i]);
                }
        }

        return 1;
}


static int test_projected_rayleigh_ritz(void)
{
        const LIS_REAL standard_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)4.0,
                (LIS_REAL)7.0
        };

        const LIS_REAL generalized_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)4.0,
                (LIS_REAL)12.0,
                (LIS_REAL)28.0
        };

        const LIS_REAL generalized_b[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)3.0,
                (LIS_REAL)4.0
        };

        LIS_MATRIX A_standard;
        LIS_MATRIX A_generalized;
        LIS_MATRIX B_generalized;

        LIS_INT err;

        A_standard = NULL;
        A_generalized = NULL;
        B_generalized = NULL;

        err = create_projector_diagonal(
                standard_a,
                &A_standard);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_a,
                &A_generalized);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_b,
                &B_generalized);

        if( err ) goto fail;

        if( run_rayleigh_ritz_case(
                    "standard",
                    A_standard,
                    NULL) )
        {
                goto fail;
        }

        if( run_rayleigh_ritz_case(
                    "generalized",
                    A_generalized,
                    B_generalized) )
        {
                goto fail;
        }

        lis_matrix_destroy(
                B_generalized);

        lis_matrix_destroy(
                A_generalized);

        lis_matrix_destroy(
                A_standard);

        printf(
                "PASS: native FEAST projected Rayleigh-Ritz "
                "standard/generalized\n");

        return 0;

fail:

        if( B_generalized!=NULL )
        {
                lis_matrix_destroy(
                        B_generalized);
        }

        if( A_generalized!=NULL )
        {
                lis_matrix_destroy(
                        A_generalized);
        }

        if( A_standard!=NULL )
        {
                lis_matrix_destroy(
                        A_standard);
        }

        return 1;
}


#define CERT_TEST_RANK 4


static int set_certification_unit_vector(
        LIS_VECTOR v,
        LIS_INT index)
{
        LIS_INT err;

        err = lis_vector_set_all(
                0.0,
                v);

        if( err )
                return err;

        err = lis_vector_set_value(
                LIS_INS_VALUE,
                index,
                1.0,
                v);

        if( err )
                return err;

        return LIS_SUCCESS;
}


static int check_certified_interval_pair(
        const char *name,
        LIS_SCALAR value[],
        LIS_INT selected_index[],
        LIS_INT selected_count,
        LIS_REAL residual[],
        LIS_REAL trace,
        LIS_REAL trace_error,
        LIS_REAL max_residual,
        LIS_INT converged,
        LIS_INT expect_converged)
{
        LIS_REAL first;
        LIS_REAL second;
        LIS_REAL low;
        LIS_REAL high;

        if( selected_count!=2 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s selected_count=%d expected 2\n",
                        name,
                        (int)selected_count);

                return 1;
        }

        if( selected_index[0]<0 ||
            selected_index[0]>=CERT_TEST_RANK ||
            selected_index[1]<0 ||
            selected_index[1]>=CERT_TEST_RANK )
        {
                fprintf(
                        stderr,
                        "FAIL: %s invalid selected index\n",
                        name);

                return 1;
        }

        first =
                projector_scalar_real(
                        value[
                                selected_index[0]]);

        second =
                projector_scalar_real(
                        value[
                                selected_index[1]]);

        if( first<second )
        {
                low = first;
                high = second;
        }
        else
        {
                low = second;
                high = first;
        }

        if( fabs((double)(low-2.0))>1.0e-12 ||
            fabs((double)(high-4.0))>1.0e-12 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s selected values %.17g %.17g "
                        "expected 2 and 4\n",
                        name,
                        (double)low,
                        (double)high);

                return 1;
        }

        if( residual[0]>1.0e-12 ||
            residual[1]>1.0e-12 ||
            max_residual>1.0e-12 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s residual certification "
                        "%.17g %.17g max %.17g\n",
                        name,
                        (double)residual[0],
                        (double)residual[1],
                        (double)max_residual);

                return 1;
        }

        if( fabs((double)(trace-6.0))>1.0e-12 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s trace %.17g expected 6\n",
                        name,
                        (double)trace);

                return 1;
        }

        if( !isfinite((double)trace_error) )
        {
                fprintf(
                        stderr,
                        "FAIL: %s non-finite trace error\n",
                        name);

                return 1;
        }

        if( converged!=expect_converged )
        {
                fprintf(
                        stderr,
                        "FAIL: %s converged=%d expected %d\n",
                        name,
                        (int)converged,
                        (int)expect_converged);

                return 1;
        }

        return 0;
}


static int run_ritz_certification_case(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B)
{
        LIS_VECTOR vectors[CERT_TEST_RANK];

        LIS_SCALAR values[CERT_TEST_RANK];

#ifdef _COMPLEX
        LIS_SCALAR bad_values[CERT_TEST_RANK];
#endif

        LIS_INT selected_index[CERT_TEST_RANK];
        LIS_REAL residual[CERT_TEST_RANK];

        LIS_INT selected_count;
        LIS_INT converged;

        LIS_REAL trace;
        LIS_REAL trace_error;
        LIS_REAL max_residual;

        LIS_INT i;
        LIS_INT err;

        for(i=0;i<CERT_TEST_RANK;i++)
        {
                vectors[i] = NULL;
                selected_index[i] = -1;
                residual[i] = 0.0;
        }

        values[0] = (LIS_SCALAR)1.0;
        values[1] = (LIS_SCALAR)2.0;
        values[2] = (LIS_SCALAR)4.0;
        values[3] = (LIS_SCALAR)7.0;

        for(i=0;i<CERT_TEST_RANK;i++)
        {
                err = lis_vector_duplicate(
                        A,
                        &vectors[i]);

                if( err ) goto fail;

                err = set_certification_unit_vector(
                        vectors[i],
                        i);

                if( err ) goto fail;
        }

        /*
         * First pass:
         * interval and physical residuals are valid, but there is
         * no previous trace, so FEAST must not report convergence.
         */
        selected_count = -1;
        converged = LIS_TRUE;
        trace = -1.0;
        trace_error = -1.0;
        max_residual = -1.0;

        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                -1,
                0.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s first certification returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        if( check_certified_interval_pair(
                    name,
                    values,
                    selected_index,
                    selected_count,
                    residual,
                    trace,
                    trace_error,
                    max_residual,
                    converged,
                    LIS_FALSE) )
        {
                goto fail;
        }

        /*
         * Second identical pass:
         * stable count + stable trace + certified residuals.
         */
        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                2,
                (LIS_REAL)6.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s second certification returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        if( check_certified_interval_pair(
                    name,
                    values,
                    selected_index,
                    selected_count,
                    residual,
                    trace,
                    trace_error,
                    max_residual,
                    converged,
                    LIS_TRUE) )
        {
                goto fail;
        }

        if( trace_error>1.0e-12 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s stable trace error %.17g\n",
                        name,
                        (double)trace_error);

                goto fail;
        }

        /*
         * A trace change must keep the outer FEAST iteration open.
         */
        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                2,
                (LIS_REAL)5.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err ) goto fail;

        if( converged!=LIS_FALSE ||
            trace_error<=1.0e-12 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s trace-change convergence gate\n",
                        name);

                goto fail;
        }

        /*
         * A count change must also keep the iteration open even
         * when the trace itself is unchanged.
         */
        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                1,
                (LIS_REAL)6.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err ) goto fail;

        if( converged!=LIS_FALSE )
        {
                fprintf(
                        stderr,
                        "FAIL: %s count-change convergence gate\n",
                        name);

                goto fail;
        }

        /*
         * Closed interval convention:
         * lambda=2 and lambda=4 remain selected at the boundaries.
         */
        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)2.0,
                (LIS_REAL)4.0,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                2,
                (LIS_REAL)6.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err ) goto fail;

        if( check_certified_interval_pair(
                    name,
                    values,
                    selected_index,
                    selected_count,
                    residual,
                    trace,
                    trace_error,
                    max_residual,
                    converged,
                    LIS_TRUE) )
        {
                fprintf(
                        stderr,
                        "FAIL: %s closed interval contract\n",
                        name);

                goto fail;
        }

        /*
         * Empty interval is a valid FEAST result.  The first empty
         * pass cannot converge, but a second stable empty pass can.
         */
        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)4.5,
                (LIS_REAL)6.0,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                -1,
                0.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err ) goto fail;

        if( selected_count!=0 ||
            trace!=0.0 ||
            max_residual!=0.0 ||
            converged!=LIS_FALSE )
        {
                fprintf(
                        stderr,
                        "FAIL: %s first empty interval pass\n",
                        name);

                goto fail;
        }

        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)4.5,
                (LIS_REAL)6.0,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                values,
                vectors,
                0,
                0.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err ) goto fail;

        if( selected_count!=0 ||
            trace!=0.0 ||
            trace_error!=0.0 ||
            max_residual!=0.0 ||
            converged!=LIS_TRUE )
        {
                fprintf(
                        stderr,
                        "FAIL: %s stable empty interval pass\n",
                        name);

                goto fail;
        }

#ifdef _COMPLEX

        /*
         * Hermitian FEAST must not silently accept a materially
         * non-real Ritz value.
         */
        for(i=0;i<CERT_TEST_RANK;i++)
        {
                bad_values[i] =
                        values[i];
        }

        bad_values[1] +=
                (LIS_SCALAR)(
                        1.0e-6
                        *
                        _Complex_I);

        err = lis_efeast_certify_ritz_pairs(
                A,
                B,
                (LIS_REAL)1.5,
                (LIS_REAL)4.5,
                (LIS_REAL)1.0e-12,
                CERT_TEST_RANK,
                bad_values,
                vectors,
                -1,
                0.0,
                selected_index,
                residual,
                &selected_count,
                &trace,
                &trace_error,
                &max_residual,
                &converged);

        if( err!=LIS_BREAKDOWN )
        {
                fprintf(
                        stderr,
                        "FAIL: %s complex non-real Ritz gate "
                        "returned %d expected %d\n",
                        name,
                        (int)err,
                        (int)LIS_BREAKDOWN);

                goto fail;
        }

#endif

        printf(
                "PASS: %s FEAST interval/residual/trace certification\n",
                name);

        for(i=0;i<CERT_TEST_RANK;i++)
        {
                if( vectors[i]!=NULL )
                {
                        lis_vector_destroy(
                                vectors[i]);
                }
        }

        return 0;

fail:

        for(i=0;i<CERT_TEST_RANK;i++)
        {
                if( vectors[i]!=NULL )
                {
                        lis_vector_destroy(
                                vectors[i]);
                }
        }

        return 1;
}


static int test_ritz_certification(void)
{
        const LIS_REAL standard_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)4.0,
                (LIS_REAL)7.0
        };

        const LIS_REAL generalized_a[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)4.0,
                (LIS_REAL)12.0,
                (LIS_REAL)28.0
        };

        const LIS_REAL generalized_b[PROJECTOR_TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)2.0,
                (LIS_REAL)3.0,
                (LIS_REAL)4.0
        };

        LIS_MATRIX A_standard;
        LIS_MATRIX A_generalized;
        LIS_MATRIX B_generalized;

        LIS_INT err;

        A_standard = NULL;
        A_generalized = NULL;
        B_generalized = NULL;

        err = create_projector_diagonal(
                standard_a,
                &A_standard);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_a,
                &A_generalized);

        if( err ) goto fail;

        err = create_projector_diagonal(
                generalized_b,
                &B_generalized);

        if( err ) goto fail;

        if( run_ritz_certification_case(
                    "standard",
                    A_standard,
                    NULL) )
        {
                goto fail;
        }

        if( run_ritz_certification_case(
                    "generalized",
                    A_generalized,
                    B_generalized) )
        {
                goto fail;
        }

        lis_matrix_destroy(
                B_generalized);

        lis_matrix_destroy(
                A_generalized);

        lis_matrix_destroy(
                A_standard);

        printf(
                "PASS: native FEAST interval/residual/trace "
                "certification standard/generalized\n");

        return 0;

fail:

        if( B_generalized!=NULL )
        {
                lis_matrix_destroy(
                        B_generalized);
        }

        if( A_generalized!=NULL )
        {
                lis_matrix_destroy(
                        A_generalized);
        }

        if( A_standard!=NULL )
        {
                lis_matrix_destroy(
                        A_standard);
        }

        return 1;
}


int main(
        int argc,
        char **argv)
{
        LIS_ESOLVER esolver;
        LIS_INT err;

        esolver = NULL;

        err = lis_initialize(
                &argc,
                &argv);

        if( err )
        {
                return (int)err;
        }

        if( test_gauss_legendre() )
        {
                goto fail;
        }

        if( test_feast_reference_contour() )
        {
                goto fail;
        }

        if( test_contour_symmetry() )
        {
                goto fail;
        }

        if( test_spectral_projector() )
        {
                goto fail;
        }

        if( test_filtered_block_basis() )
        {
                goto fail;
        }

        if( test_projected_rayleigh_ritz() )
        {
                goto fail;
        }

        if( test_ritz_certification() )
        {
                goto fail;
        }

        err = lis_esolver_create(
                &esolver);

        if( err )
        {
                goto fail;
        }

        err = lis_esolver_set_option(
                "-e feast "
                "-ss 12 "
                "-feast_emin 1.0 "
                "-feast_emax 5.0 "
                "-feast_ratio 0.30 "
                "-feast_points 8",
                esolver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: valid contour options rejected\n");

                goto fail;
        }

        err = lis_efeast_check_params(
                esolver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: valid contour configuration rejected\n");

                goto fail;
        }

        lis_esolver_destroy(esolver);
        lis_finalize();

        printf(
                "PASS: native FEAST Gauss contour\n");

        return 0;

fail:

        if( esolver!=NULL )
        {
                lis_esolver_destroy(esolver);
        }

        lis_finalize();

        return 1;
}
