# LIS BUCKLE – etap 70 / RC1 prototypu (kandydat poprawki 59)

Stan: **samodzielny kandydat do integracji, NIE wydanie produkcyjne LIS/CCX.**
Kod C11 etapów 55–69, niezależny certyfikator residuum oryginalnego problemu własnego, pełna regresja. Bez BLAS, OpenBLAS, LAPACK, METIS i SuiteSparse. Wersja wykorzystuje libm, libc i libgomp. Adapter domyślnie używa jednego wątku; TSan pozostaje nierozstrzygnięty (log etapu 68 zachowany, nie uruchomiono nowego TSan).

## Nowy interfejs
`mode_certificate.h` / `mode_certificate.c`: p70_certify oblicza
`||Kx-lambda Gx||_2 / (||Kx||_2 + |lambda| ||Gx||_2)`
na oryginalnych CSR, niezależnie od faktorów i przesunięcia. Używa akumulacji long double i normalizuje x według maksymalnego modułu. Akceptacja wg podanego progu, np. 1e-10. Odrzuca niepoprawne dane, wektor zerowy, nieokreślony iloraz 0/0 i przekroczenie budżetu dwóch tablic long double. P55_OK oznacza poprawne obliczenie, nie automatycznie akceptację; sprawdź `accepted`. Dane wejściowe i długości faktycznych tablic muszą być poprawne po stronie wywołującego. Współdzieli walidację CSR z adapterem etapu 69, nie korzysta z przesuniętego operatora ani faktorów. Certyfikat residuum sam w sobie nie dowodzi kompletności zbioru modów ani poprawności ich krotności. Nie jest to identyczna definicja residuum z etapów 52/54, jeśli używały innego mianownika.

`test_certificate.c`: 864 przypadki analityczne (432 poprawne pary i 432 celowo błędne wartości własne), 12 odzyskanych par (dwa testy shift-invert po 6 modów) oraz 8 błędnych danych. Problem 12 DOF, G dodatnio określona, wartości własne 1..12, obrócone pary współrzędnych, skale 1e-100..1e100. Przesunięcia 0.25 i 5.30. Iteracja testowa z deflacją w iloczynie skalarnym G służy wyłącznie do weryfikacji tego problemu; nie zastępuje eigensolvera LIS ani nie jest ogólnym algorytmem dla nieokreślonej G. Wyniki: 884/884 PASS w obu konfiguracjach, max residual 1.457e-16 dla poprawnych par, max eigenvalue error 8.882e-16.

## Budowanie i regresja
```
make -j4 stage70-check
```
GCC 14.2, Clang 17, Linux x86_64. Dla sanitizerów z OpenMP Makefile używa GCC. Pełne `make -j4 stage69-check` zakończyło się kodem 0 (logs/regression.txt, logs/regression.exit). Po wyodrębnieniu walidacji CSR adapter ponownie przetestowano w obu konfiguracjach (logs/bridge-final.txt). `make certificate` i `make certificate-sanitize` zakończone poprawnie (osobne logi). Wszystkie 10177/10177 standardowych, 10175/10175 sanitizowanych oraz 162/162 portable testów PASS, razem 20514 wykonań; powtórzenie bridge nie doliczone. Różna liczba testów parallel (526 vs 524) wynika z dynamicznie ustalanych progów budżetu. Czasy w logach regresji równoległej kompilacji nie są miarodajnym benchmarkiem. Nie przeprowadzono TSan w tym etapie, MPI, ABI ani make check oryginalnego LIS. Wcześniejsze logi i wyniki dotyczą wcześniejszych wersji, nie tej regresji.

## Integracja i warunki wydania
1. Dostarczyć zmodyfikowane źródła LIS etapu 52/54, wersję i źródła CCX, Makefile/konfigurację, model .inp i wszystkie INCLUDE. Alternatywnie eksport K/G, RHS i modów może umożliwić oddzielną walidację macierzową, ale nie zastąpi pełnego BUCKLE.
2. Zweryfikować konwencję znaku, sigma, typy LIS_INT/LIS_SCALAR, zakres i własność CSR; nie rzutować niezgodnych typów wskaźników. Nie interpretować bloku lokalnego MPI jako macierzy globalnej. Użyć p69_factorize / p56_apply zgodnie z README-stage69.md. Nie używać zachowanego starego faktora po błędzie dla nowego operatora.
3. Domyślnie jeden wątek; OpenMP tylko eksperymentalnie, do czasu rozwiązania ostrzeżeń TSan. Automatyczny cache wyłączony w adapterze; w razie integracji klucz musi obejmować operator, sigma i konfigurację.
4. Uruchomić identyczny model z bazowym etapem 52 i kandydatem, minimum 5 powtórzeń w tych samych warunkach (osobno cold/warm cache, ten sam sprzęt i liczba wątków). Zmierzyć wall time, szczytowy RSS całego CCX, osobno eliminację i preconditioner. Nie sumować szczytów różnych faz.
5. Sprawdzić wszystkie sześć modów w oryginalnym problemie, wartości własne oraz odpowiadające przestrzenie własne dla modów wielokrotnych; porównać residuum wg tej samej definicji co baseline. Zbadać tolerancje i przypadki osobliwe.
6. Przebudować LIS/CCX od zera i sprawdzić make check, MPI i ABI; przeprowadzić weryfikację wyścigów z kompatybilnym runtime OpenMP. Dopiero po tym decyzja produkcyjna.

## Status
Brak oryginalnych źródeł LIS 52/54 i rzeczywistego modelu CCX 152592 DOF w dostępnych danych. Nie zmodyfikowano LIS ani CCX, nie opublikowano na GitHubie. Nie potwierdzono celu 3–5 s BUCKLE / <500 MiB. Etap 52 pozostaje referencyjną wersją na podstawie wcześniejszych wyników przekazanych przez użytkownika. Nie planuje się dalszego zwiększania liczby etapów samodzielnego prototypu bez danych integracyjnych.
