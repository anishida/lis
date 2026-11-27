# LIS BUCKLE – etap 63 / kandydat poprawki 52

Samodzielny prototyp C11 oparty na etapie 62. Nie jest poprawką do oficjalnego LIS. Bez BLAS, LAPACK, METIS, SuiteSparse. Wszystkie moduły poprzednich etapów pozostają dostępne.

## API

`minimum_degree.h`: `p63_analyze` wyznacza dokładne Minimum Degree dla dolnego CSC (wiersze rosnące, bez duplikatów, jawne zera traktowane strukturalnie). `perm[k]` = oryginalny indeks zmiennej na pozycji k; `inv` = odwrotność. `p63_free` zwalnia wynik. Ustaw `*out = NULL` przed pierwszym wywołaniem. Wynik zachowany na błędzie. Remisy deterministyczne według najmniejszego pierwotnego indeksu. Graf jawny z dynamicznymi listami sąsiedztwa i indeksowanym kopcem minimum; po usunięciu wierzchołka sąsiedzi tworzą klikę. To exact MD, nie AMD ani algorytm grafu ilorazowego.

`p63_factorize` buduje uporządkowane dolne CSC, wywołuje p62/p61 i składa permutację porządkowania z permutacją pivotów numerycznych. `p56_apply` przyjmuje i zwraca wektory w pierwotnej numeracji. Statystyki oraz poprzedni kontekst są zachowywane na błędzie. Pole `summary` w statystykach zawiera wyłącznie liczby, wskaźniki perm/inv są NULL. Przykładowe wywołania w test_ordering.c.

## Limity

`max_work_bytes` ogranicza jawne alokacje etapu porządkowania (w tym jednoczesne stare/nowe bufory przy rozszerzaniu) oraz oddzielnie transformacji CSC. `max_slots` ogranicza przewidywaną liczbę miejsc poza diagonalą L; `max_pair_checks` ogranicza sprawdzenia par w klikach, nie wszystkie operacje grafowe. Budżet nie obejmuje metadanych alokatora, ewentualnego wewnętrznego bufora qsort, starego kontekstu, danych wejściowych i sumy równoczesnych alokacji modułów p62/p61. Nie jest gwarancją RSS. Symboliczny rozmiar frontu oraz liczba miejsc w L mogą odbiegać od numerycznych z powodu pivotów opóźnionych i zer numerycznych. Graf jawny może być kosztowny przy wielkich klikach. MD nie gwarantuje zmniejszenia wypełnienia ani maksymalnego frontu dla każdej macierzy. Brak automatycznego wyboru kolejności lub fallbacku do naturalnej.

## Testy i pomiary

`make ordering` – 903 testy; `make ordering-sanitize` – ASan/UBSan; `make -j4 stage63-check` – regresja 5208 przypadków na konfigurację, obie konfiguracje przechodzą. 700 porównań z niezależnym gęstym grafem (dokładna kolejność, stopnie, liczba miejsc, fronty, liczba par), w tym 300 rozwiązań; 160 układów nieokreślonych z pivotami 2x2; 40 przypadków walidacji/limitów; 3 duże układy. Maksymalny błąd rozwiązania 1.554e-15, komponentowy błąd wsteczny 7.901e-16. Surowe logi w test-results.txt, także oddzielny przebieg zoptymalizowany.

| Siatka | naturalne miejsca L | MD miejsca L | front naturalny | front MD |
|---|---:|---:|---:|---:|
| 40x40 | 62439 | 19904 | 41 | 63 |
| 128x128 | 2080895 | 400929 | 129 | 244 |

Pomiar CPU w pojedynczym zoptymalizowanym przebiegu dla 40x40: naturalna analiza i faktoryzacja 0.557915 s, MD wraz z porządkowaniem 0.107732 s. Pamięć robocza modułu numerycznego 1707248 -> 707188 bajtów. Nie jest to RSS. Dla 128x128 wykonano tylko analizę symboliczną (MD 0.102786 s). Dla trójdiagonalnego 152592 DOF nie ma oszczędności miejsc w L, naturalna 0.136163 s, MD 0.150791 s. Czasy pojedyncze, zależne od środowiska, nie pomiary CCX.

## Następne kroki

Porównanie kolejności i wybór na podstawie przewidywanych miejsc L i maksymalnego frontu; redukcja kosztu jawnego grafu; uporządkowanie separatorowe dla siatek; numeryczne łączenie superwęzłów; pełna integracja LIS. Brakuje oryginalnego kodu etapów 52 i 54 oraz macierzy CCX, nie wykonano testów pełnego BUCKLE, MPI, ABI. Nie publikowano zmian na GitHubie.
