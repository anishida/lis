#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "lis.h"
#include "lis_esolver.h"


static int expect_failure(
        LIS_INT err,
        const char *label)
{
        if( err==LIS_SUCCESS )
        {
                fprintf(
                        stderr,
                        "FAIL: %s unexpectedly succeeded\n",
                        label);
                return 1;
        }

        return 0;
}


int main(
        int argc,
        char **argv)
{
        LIS_ESOLVER esolver;
        LIS_INT err;
        LIS_INT nesolver;
        char name[128];

        esolver = NULL;

        err = lis_initialize(&argc,&argv);
        if( err ) return (int)err;

        err = lis_esolver_create(&esolver);

        if( err )
        {
                lis_finalize();
                return (int)err;
        }

        err = lis_esolver_set_option(
                "-e feast "
                "-ss 12 "
                "-feast_emin 1.0 "
                "-feast_emax 5.0 "
                "-feast_points 8",
                esolver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: valid native FEAST options rejected: %d\n",
                        (int)err);
                goto fail;
        }

        err = lis_esolver_get_esolver(
                esolver,
                &nesolver);

        if( err ) goto fail;

        if( nesolver!=LIS_ESOLVER_FEAST )
        {
                fprintf(
                        stderr,
                        "FAIL: -e feast selected %d, expected %d\n",
                        (int)nesolver,
                        (int)LIS_ESOLVER_FEAST);
                goto fail;
        }

        err = lis_esolver_get_esolvername(
                nesolver,
                name);

        if( err ) goto fail;

        if( strcmp(name,"FEAST")!=0 )
        {
                fprintf(
                        stderr,
                        "FAIL: FEAST name is '%s'\n",
                        name);
                goto fail;
        }

        err = lis_efeast_check_params(esolver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: valid FEAST configuration rejected: %d\n",
                        (int)err);
                goto fail;
        }

        err = lis_esolver_set_option(
                "-feast_points 0",
                esolver);

        if( expect_failure(
                err,
                "-feast_points 0") )
        {
                goto fail;
        }

        err = lis_esolver_set_option(
                "-feast_points 8",
                esolver);

        if( err ) goto fail;

        err = lis_esolver_set_option(
                "-feast_emin 5.0 "
                "-feast_emax 1.0",
                esolver);

        if( err ) goto fail;

        err = lis_efeast_check_params(esolver);

        if( expect_failure(
                err,
                "reversed FEAST interval") )
        {
                goto fail;
        }

        err = lis_esolver_set_option(
                "-feast_emin 1.0 "
                "-feast_emax 5.0",
                esolver);

        if( err ) goto fail;

        err = lis_efeast_check_params(esolver);

        if( err ) goto fail;

        lis_esolver_destroy(esolver);
        lis_finalize();

        printf("PASS: native FEAST foundation\n");

        return 0;

fail:

        if( esolver!=NULL )
        {
                lis_esolver_destroy(esolver);
        }

        lis_finalize();

        return 1;
}
