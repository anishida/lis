# LIS BUCKLE – etap 68, kandydat poprawki 57

Samodzielny prototyp C11, kontynuacja etapu 67. Nie jest poprawką zintegrowaną z LIS/CCX. Bez BLAS, OpenBLAS, LAPACK, METIS i SuiteSparse. Równoległa ścieżka wymaga OpenMP.

## Zmiany
- `parallel_ldlt.c`: rozpoznanie łańcucha w zakontraktowanym drzewie przed alokacją harmonogramu, przejście do p65_numeric z zachowaniem walidacji. Odroczone alokacje buforów numerycznych i tablic składania, wcześniejsze zwalnianie buforów. Budżet per wywołanie chroniony atomic_flag C11, alokator poza sekcją krytyczną, status błędu atomic_int. Zachowano API p67.
- `factor_cache.c`: zbiorczy zapis i odczyt tablic binary64 w natywnym little endian; zbiorczy zapis indeksów i permutacji. Zachowano format LDLTC001, CRC32 i walidację. P68_CACHE_SCALAR wymusza poprzednią ścieżkę element po elemencie; na big endian pozostaje konwersja. Bez zmiany semantyki rename/fsync; katalog nie jest synchronizowany.
- `test_memory_cache.c`: test dokładnego budżetu łańcucha (i budżetu o bajt mniejszego), porównanie bitowe rozwiązań i plików, odczyt krzyżowy obu wariantów cache, 60 współbieżnych niezależnych faktoryzacji z zagnieżdżonym OpenMP. Nie zmieniać współdzielonego trybu SIMD podczas faktoryzacji; wymiana tego samego kontekstu wymaga synchronizacji wywołującego.

## Kompilacja i testy
```
make -j2 stage68-check
make parallel-bench
make memory-cache
```
GCC 14.2 (OpenMP, również ASan/UBSan dla testów równoległych), Clang 17 (pozostałe sanitizery). W tym przebiegu stage67-check zakończony kodem 0, memory-cache i memory-cache-sanitize zakończone kodem 0. 8857 przypadków podstawowych w każdej konfiguracji (8259 z etapów 55–66, 526 OpenMP, 72 nowe), plus osobny przebieg 162 przypadków bez SIMD. Liczba przypadków OpenMP zależy od pętli limitów pamięci, nie jest stała. Dwa zestawy testów zapisują czasy, lecz benchmark-parallel.txt i benchmark-cache.txt są niezależnymi pomiarami wykonanymi kolejno po zakończeniu regresji.

## Pamięć i wydajność
Syntetyczna macierz trójdiagonalna n=152592, żądane 4 wątki: 41047256 bajtów w etapie 67, 23347608 bajtów w etapie 68 (43,12% mniej). Dokładnie tyle samo co przy jednym wątku w benchmarku etapu 68. Dla siatki 128x128: 13720852 -> 12785408 bajtów; bloki 8192 DOF: 8539304 -> 8047424 bajtów. Pamięć raportowana to monitorowane alokacje modułu, nie całe RSS. Wyniki etapu 67 pochodzą z jego wcześniejszego benchmarku; pomiary nie są jednoczesne. Budżet nie obejmuje nagłówków alokatora, stosów, runtime OpenMP, macierzy wejściowej, analizy symbolicznej i końcowego/starego adaptera. Rezerwacje pamięci liczone przed calloc, więc szczyt może obejmować oczekującą alokację.

Cache: w dwóch testach pliki obu trybów są identyczne bajtowo, również z permutacją i pivotami 2x2. Mediany 3 powtórzeń na 17313852 bajtach: zapis 0.059212 -> 0.056227 s, odczyt 0.052187 -> 0.046984 s. Na 1120316 bajtach: zapis 0.004655 -> 0.004285 s, odczyt 0.003289 -> 0.002890 s. Są to pomiary lokalnego cache plików, nie gwarancja wydajności dysku ani CCX; fsync w zapisie pozostaje. Nie testowano big endian, tylko wymuszoną ścieżkę skalarną na little endian. Klucz 32 bajty nadal dostarcza wywołujący.

## Bezpieczeństwo wielowątkowe – blokada wydania
ASan/UBSan i testy bitowe przeszły, ale ThreadSanitizer GCC/libgomp nadal zgłasza konflikt między zapisem wyniku frontu a odczytem po niejawnej barierze omp for (parallel_ldlt.c:214/217). Uruchomiono halt_on_error=1, kod 66. Nie rozstrzygnięto, czy to ograniczenie instrumentacji libgomp czy rzeczywisty wyścig; brak libomp/Archer w środowisku. Nie wyciszano ostrzeżenia. Należy zweryfikować na runtime z obsługą TSan/Archer i przeprowadzić audyt synchronizacji przed wydaniem równoległym. Zmiana blokady budżetu nie dowodzi braku wyścigów w pozostałym kodzie.

## Kolejne etapy
69: integracja ze zmodyfikowanym LIS/CCX; 70: pełna regresja, benchmark i wydanie. Potrzebne oryginalne źródła etapów 52/54 (albo odpowiadające im patche oraz baza) i rzeczywisty przypadek CCX 152592 DOF. Nie wykonano benchmarku rzeczywistego CCX, pełnego make check LIS, MPI, ABI ani publikacji GitHub. Nie traktować jako wersji produkcyjnej. Nie obiecuje się osiągnięcia docelowego czasu BUCKLE 3–5 s ani RAM poniżej 500 MiB.
