CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Werror
SAN_CC ?= clang
.PHONY: check sanitize clean
check: test_packed
	./test_packed
test_packed: factor_private.h packed_ldlt.c packed_ldlt.h test_packed.c
	$(CC) $(CFLAGS) packed_ldlt.c test_packed.c -lm -o $@
sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer packed_ldlt.c test_packed.c -lm -o test_packed_san
	./test_packed_san
clean:
	rm -f test_qmd_plan test_qmd_plan_san test_schur_omp test_schur_omp_san test_certificate test_certificate_san test_bridge test_bridge_san test_memory_cache test_memory_cache_san cache_scalar.o cache_scalar_san.o test68_tsan test_parallel test_parallel_san test_simd test_simd_san test_simd_portable test_supernodal test_supernodal_san test_nd test_nd_san test_ordering test_ordering_san test_symbolic test_symbolic_san test_multifrontal test_multifrontal_san test_delayed test_delayed_san test_bk test_bk_san test_pivot test_pivot_san test_cache test_cache_san test_packed test_packed_san test_adapter test_adapter_san
.PHONY: adapter adapter-sanitize all-check
adapter: test_adapter
	./test_adapter
test_adapter: factor_private.h packed_ldlt.c packed_ldlt.h lis_adapter.c lis_adapter.h test_adapter.c
	$(CC) $(CFLAGS) packed_ldlt.c lis_adapter.c test_adapter.c -lm -o $@
adapter-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer packed_ldlt.c lis_adapter.c test_adapter.c -lm -o test_adapter_san
	./test_adapter_san
all-check: check adapter sanitize adapter-sanitize
.PHONY: cache cache-sanitize full-check
cache: test_cache
	./test_cache
test_cache: factor_cache.c factor_cache.h factor_private.h packed_ldlt.c lis_adapter.c test_cache.c
	$(CC) $(CFLAGS) packed_ldlt.c lis_adapter.c factor_cache.c test_cache.c -lm -o $@
cache-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer packed_ldlt.c lis_adapter.c factor_cache.c test_cache.c -lm -o test_cache_san
	./test_cache_san
full-check: all-check cache cache-sanitize
.PHONY: pivot pivot-sanitize stage58-check
pivot: test_pivot
	./test_pivot
test_pivot: pivot_ldlt.c pivot_ldlt.h packed_ldlt.c lis_adapter.c factor_cache.c test_pivot.c
	$(CC) $(CFLAGS) packed_ldlt.c lis_adapter.c factor_cache.c pivot_ldlt.c test_pivot.c -lm -o $@
pivot-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer packed_ldlt.c lis_adapter.c factor_cache.c pivot_ldlt.c test_pivot.c -lm -o test_pivot_san
	./test_pivot_san
stage58-check: full-check pivot pivot-sanitize
.PHONY: bk bk-sanitize stage59-check
bk: test_bk
	./test_bk
test_bk: bk_ldlt.c bk_ldlt.h pivot_ldlt.c packed_ldlt.c lis_adapter.c factor_cache.c test_bk.c
	$(CC) $(CFLAGS) packed_ldlt.c lis_adapter.c factor_cache.c pivot_ldlt.c bk_ldlt.c test_bk.c -lm -o $@
bk-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer packed_ldlt.c lis_adapter.c factor_cache.c pivot_ldlt.c bk_ldlt.c test_bk.c -lm -o test_bk_san
	./test_bk_san
stage59-check: stage58-check bk bk-sanitize
.PHONY: delayed delayed-sanitize stage60-check
DELAY_SRC=packed_ldlt.c lis_adapter.c factor_cache.c pivot_ldlt.c bk_ldlt.c delayed_ldlt.c test_delayed.c
delayed: test_delayed
	./test_delayed
test_delayed: $(DELAY_SRC) delayed_ldlt.h
	$(CC) $(CFLAGS) $(DELAY_SRC) -lm -o $@
delayed-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(DELAY_SRC) -lm -o test_delayed_san
	./test_delayed_san
stage60-check: stage59-check delayed delayed-sanitize
.PHONY: multifrontal multifrontal-sanitize stage61-check
MF_SRC=packed_ldlt.c lis_adapter.c factor_cache.c multifrontal_ldlt.c schur_simd.c test_multifrontal.c
multifrontal: test_multifrontal
	./test_multifrontal
test_multifrontal: $(MF_SRC) multifrontal_ldlt.h
	$(CC) $(CFLAGS) $(MF_SRC) -lm -o $@
multifrontal-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(MF_SRC) -lm -o test_multifrontal_san
	./test_multifrontal_san
stage61-check: stage60-check multifrontal multifrontal-sanitize
.PHONY: symbolic symbolic-sanitize stage62-check
SYM_SRC=packed_ldlt.c lis_adapter.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c test_symbolic.c
symbolic: test_symbolic
	./test_symbolic
test_symbolic: $(SYM_SRC) symbolic_ldlt.h
	$(CC) $(CFLAGS) $(SYM_SRC) -lm -o $@
symbolic-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(SYM_SRC) -lm -o test_symbolic_san
	./test_symbolic_san
stage62-check: stage61-check symbolic symbolic-sanitize
.PHONY: ordering ordering-sanitize stage63-check
ORD_SRC=packed_ldlt.c lis_adapter.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c minimum_degree.c test_ordering.c
ordering: test_ordering
	./test_ordering
test_ordering: $(ORD_SRC) minimum_degree.h
	$(CC) $(CFLAGS) $(ORD_SRC) -lm -o $@
ordering-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(ORD_SRC) -lm -o test_ordering_san
	./test_ordering_san
stage63-check: stage62-check ordering ordering-sanitize
.PHONY: nd nd-sanitize stage64-check
ND_SRC=packed_ldlt.c lis_adapter.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c minimum_degree.c nested_dissection.c test_nd.c
nd: test_nd
	./test_nd
test_nd: $(ND_SRC) nested_dissection.h
	$(CC) $(CFLAGS) $(ND_SRC) -lm -o $@
nd-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(ND_SRC) -lm -o test_nd_san
	./test_nd_san
stage64-check: stage63-check nd nd-sanitize
.PHONY: supernodal supernodal-sanitize stage65-check
SUPER_SRC=packed_ldlt.c lis_adapter.c factor_cache.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c minimum_degree.c nested_dissection.c test_supernodal.c
supernodal: test_supernodal
	./test_supernodal
test_supernodal: $(SUPER_SRC) supernodal_ldlt.h
	$(CC) $(CFLAGS) $(SUPER_SRC) -lm -o $@
supernodal-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(SUPER_SRC) -lm -o test_supernodal_san
	./test_supernodal_san
stage65-check: stage64-check supernodal supernodal-sanitize
.PHONY: simd simd-sanitize simd-portable simd-bench stage66-check
SIMD_SRC=packed_ldlt.c lis_adapter.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c minimum_degree.c nested_dissection.c test_simd.c
test_simd: $(SIMD_SRC) schur_simd.h
	$(CC) $(CFLAGS) $(SIMD_SRC) -lm -o $@
simd: test_simd
	./test_simd
simd-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(SIMD_SRC) -lm -o test_simd_san
	./test_simd_san
simd-portable:
	$(CC) $(CFLAGS) -DP66_DISABLE_SIMD $(SIMD_SRC) -lm -o test_simd_portable
	./test_simd_portable
simd-bench: test_simd
	./test_simd bench
stage66-check: stage65-check simd simd-sanitize simd-portable

test_multifrontal test_symbolic test_ordering test_nd test_supernodal: schur_simd.h
.PHONY: parallel parallel-sanitize parallel-bench stage67-check
PAR_SRC=packed_ldlt.c lis_adapter.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c minimum_degree.c nested_dissection.c parallel_ldlt.c test_parallel.c
test_parallel: $(PAR_SRC) parallel_ldlt.h
	$(CC) $(CFLAGS) -fopenmp $(PAR_SRC) -lm -o $@
parallel: test_parallel
	./test_parallel
parallel-sanitize:
	$(CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer -no-pie $(PAR_SRC) -lm -o test_parallel_san
	./test_parallel_san
parallel-bench: test_parallel
	./test_parallel bench
stage67-check: stage66-check parallel parallel-sanitize
.PHONY: memory-cache memory-cache-sanitize stage68-check
P68_SRC=$(filter-out test_parallel.c,$(PAR_SRC)) factor_cache.c test_memory_cache.c
P68_SCALAR=-DP68_CACHE_SCALAR -Dp57_save=scalar_save -Dp57_load=scalar_load
memory-cache:
	$(CC) $(CFLAGS) $(P68_SCALAR) -c factor_cache.c -o cache_scalar.o
	$(CC) $(CFLAGS) -fopenmp $(P68_SRC) cache_scalar.o -lm -o test_memory_cache
	./test_memory_cache
memory-cache-sanitize:
	$(CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined $(P68_SCALAR) -c factor_cache.c -o cache_scalar_san.o
	$(CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer -no-pie $(P68_SRC) cache_scalar_san.o -lm -o test_memory_cache_san
	./test_memory_cache_san
stage68-check: stage67-check memory-cache memory-cache-sanitize
.PHONY: bridge bridge-sanitize stage69-check
BRIDGE_SRC=$(filter-out test_parallel.c,$(PAR_SRC)) buckle_bridge.c test_bridge.c
bridge: test_bridge
	./test_bridge
test_bridge: $(BRIDGE_SRC) buckle_bridge.h
	$(CC) $(CFLAGS) -fopenmp $(BRIDGE_SRC) -lm -o $@
bridge-sanitize:
	$(CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer -no-pie $(BRIDGE_SRC) -lm -o test_bridge_san
	./test_bridge_san
stage69-check: stage68-check bridge bridge-sanitize
.PHONY: certificate certificate-sanitize stage70-check
CERT_SRC=$(filter-out test_bridge.c,$(BRIDGE_SRC)) mode_certificate.c test_certificate.c
certificate:
	$(CC) $(CFLAGS) -fopenmp $(CERT_SRC) -lm -o test_certificate
	./test_certificate
certificate-sanitize:
	$(CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer -no-pie $(CERT_SRC) -lm -o test_certificate_san
	./test_certificate_san
stage70-check: stage69-check certificate certificate-sanitize

# Explicit experimental symbolic-plan and Schur tile tests. These are
# not part of the installed library or the default build.
PLAN_SRC=qmd_plan.c minimum_degree.c packed_ldlt.c lis_adapter.c multifrontal_ldlt.c schur_simd.c symbolic_ldlt.c test_qmd_plan.c
.PHONY: qmd-plan schur-omp qmd-plan-sanitize schur-omp-sanitize
qmd-plan:
	$(CC) $(CFLAGS) -DP66_DISABLE_SIMD $(PLAN_SRC) -lm -o test_qmd_plan
	./test_qmd_plan
qmd-plan-sanitize:
	$(SAN_CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -DP66_DISABLE_SIMD $(PLAN_SRC) -lm -o test_qmd_plan_san
	./test_qmd_plan_san
schur-omp:
	$(CC) $(CFLAGS) -fopenmp test_schur_omp.c schur_simd.c -lm -o test_schur_omp
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 ./test_schur_omp
schur-omp-sanitize:
	$(CC) -std=c11 -Wall -Wextra -Werror -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fopenmp -no-pie test_schur_omp.c schur_simd.c -lm -o test_schur_omp_san
	OMP_DYNAMIC=FALSE OMP_NUM_THREADS=4 ./test_schur_omp_san
