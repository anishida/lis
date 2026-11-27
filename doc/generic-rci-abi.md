# Experimental RCI compatibility proposal

Adapted from the archived ABI_POLICY.md, not an official LIS ABI freeze.
Scope: opt-in non-MPI real-double Serial/OpenMP, not compatibility with the
unavailable RC2/RC3. Public types and status constants in lis_rci.h and
lis_rci_bridge.h are experimental. The current proposal preserves all ten
symbols of the first reconstruction and uses additive extensions only.
Do not renumber statuses or methods, reorder or resize existing public
structs, change calling conventions or remove symbols in this proposal.
An ABI-breaking future change requires a distinct versioned interface or
explicit ABI change agreed with the maintainer. Layouts depend on the C/C++
ABI (size_t, long double, enum representation and alignment); no guarantee
across architectures, -fshort-enums or incompatible calling conventions.

Contexts/bridges are opaque and owned by the creator. Request arrays are
borrowed, inputs read-only, outputs must be completely written. Tokens are
context-local; stale descriptors cannot safely be written. Workers must
finish accessing buffers before cancellation or destruction. Independent
contexts have no shared mutable core state; one context is not reentrant.
A bridge borrows the assembled matrix; use matching installed LIS headers
and flags and verify LIS_RCI_CHECK_CLIENT_ABI before typed LIS calls. The
bridge copies the answer only on DONE; failed generic iterates are not
certified solutions. CG requires SPD operator and preconditioner; BiCGSTAB
and GMRES require fixed linear preconditioning, FGMRES may vary. Identity
and Jacobi are native bridge options; generic PRECON supports caller-owned
preconditioners. Convergence does not imply small forward error for an
ill-conditioned system.

Acceptance: compare installed headers, exported symbols and SONAME, run
previously compiled C/C++ clients without relinking, native regression,
nonfinite/callback/cancellation probes, and industrial diagnostics. Default
builds remain unchanged. Maintainer agreement and contributor license
confirmation remain required before a release or official ABI promise.
