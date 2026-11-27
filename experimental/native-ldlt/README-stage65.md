# LIS BUCKLE – etap 65 / kandydat poprawki 54

Samodzielny prototyp C11, kontynuacja etapu 64. Bez BLAS/LAPACK/METIS/SuiteSparse. Nie jest patchem oryginalnego LIS.

## Implementacja
- `p65_factorize` wyznacza ścisłe superwęzły za pomocą p62_analyze, dzieli je do zadanej szerokości max_width, kontraktuje łańcuchy drzewa i składa lokalne macierze z wielu kolumn CSC oraz aktualizacji dzieci.
- `p65_ordered_factorize` wykorzystuje Natural/MD/ND/AUTO oraz składa permutacje jak w etapie 64.
- `max_width=1` zachowuje dotychczasową ścieżkę p62/p64. Większe szerokości ograniczają liczbę frontów, ale nie gwarantują poprawy czasu, pamięci lub stabilności. Limity numeryczne działają jak wcześniej. Nie ma automatycznego ponowienia ze zmniejszoną szerokością.
- Bunch–Kaufman wybiera pivoty spośród pierwotnych kolumn frontu i opóźnionych kolumn dzieci. Niewyeliminowane kolumny wraz z dopełnieniem Schura przechodzą do rodzica. Transakcyjna podmiana kontekstu pozostaje zachowana.
- Pakowanie wynikowych L wciąż odbywa się według pojedynczych pivotów (szerokość 1 lub 2); format cache niezmieniony. Nie zaimplementowano blokowych jąder macierz–macierz, SIMD, równoległości ani numerycznego pakowania całego superwęzła w jeden panel.
- Symboliczne szacunki AUTO są nadal oparte na liczbie elementów kolumn L; nie modelują opóźnionych pivotów. Budżety fazowe nie obejmują całego RSS.

## Testy
`make stage65-check` – wszystkie dotychczasowe zestawy oraz test_supernodal, GCC i Clang ASan+UBSan.
`make supernodal` – nowe testy i pomiary CPU (minimum z trzech powtórzeń, rozwiązanie i weryfikacja poza mierzonym czasem). `make supernodal-sanitize` – analogiczny zestaw z sanitizerami.

783/783 nowych przypadków, 6929/6929 łącznie w każdej konfiguracji. W nowych testach max error x=1.221e-15, max componentwise backward=8.105e-16. 400 losowych macierzy w trzech skalach 1e-100, 1, 1e100: 12827 -> 4714 frontów; 150 dodatkowych testów z porządkowaniem; 200 gęstych nieokreślonych przypadków ze słabymi sprzężeniami: 3020 -> 0 opóźnień, 3020 pivotów 2x2. Testy granic budżetów, transakcji, cache, wielu dzieci jednego superwęzła oraz 9 konfiguracji benchmarku. Licznik przypadków nie jest licznikiem pojedynczych asercji ani wszystkich wywołań solvera.

Log sekwencyjny: `benchmark-serial.txt`. Log pełnej regresji: `test-results.txt` (czasy z równoległego make nie są miarodajnym benchmarkiem). `dependencies.txt` – ldd nowego binarium.

## Ograniczenia
Brak źródeł zmodyfikowanego LIS i rzeczywistej macierzy CCX 152592 DOF. Test trójdiagonalny z taką liczbą niewiadomych nie jest benchmarkiem CCX. Nie wykonano MPI/ABI, pełnego BUCKLE ani make check oryginalnego LIS. Nie publikowano kodu na GitHubie. Kolejny etap: wydajne jądra Schura i SIMD z kontrolowaną dokładnością, następnie OpenMP, pamięć i integracja.
