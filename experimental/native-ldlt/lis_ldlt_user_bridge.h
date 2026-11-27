#ifndef LIS_LDLT_USER_BRIDGE_H
#define LIS_LDLT_USER_BRIDGE_H
#include "lis.h"
/* Register once after lis_initialize, before creating/using solvers.
 * Experimental application-linked bridge; not an installed LIS API. */
LIS_INT lis_ldlt_experimental_register(void);
#endif
