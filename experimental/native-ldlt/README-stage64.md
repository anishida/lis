# LIS BUCKLE – etap 64 (kandydat poprawki 53)

Samodzielny prototyp C11 rozwijający etap 63. Nie jest poprawką do oryginalnego LIS. Brak bibliotek BLAS, LAPACK, METIS i SuiteSparse; test wykonywalny linkuje libm i libc.

## Nowy moduł

`nested_dissection.h/.c`: `p64_order` oraz `p64_factorize`. Dostępne Natural, dokładne Minimum Degree z etapu 63, heurystyczne Nested Dissection i AUTO. Dane wejściowe: dolny CSC z rosnącymi indeksami wierszy, bez duplikatów. Permutacja `perm[k]` oznacza pierwotny DOF w pozycji k. `inv` jest odwrotnością. Inicjalizuj wskaźniki wynikowe NULL; zwalniaj przez `p63_free` / `p56_destroy`. W razie błędu wskaźnik i statystyki pozostają niezmienione.

ND buduje graf CSR bez diagonali. Dwie przebieżki BFS wyznaczają warstwy; separator to wewnętrzna warstwa minimalizująca max(rozmiar lewej części, rozmiar prawej części) + 2*rozmiar separatora. Części są rozłączne bez krawędzi między nimi; separator umieszczany jest na końcu. Rozłączne składowe grupowane są na granicy najbliższej połowie. Liście do `leaf_size`, brak warstwy wewnętrznej lub `max_depth` kończą rekursję bez dodatkowej optymalizacji. Parametry testów: leaf_size=8, max_depth=64. To heurystyka separatorów BFS, nie pełna implementacja wielopoziomowego partycjonowania.

AUTO liczy dokładne strukturalne miejsca L i maksymalny front oraz proxy kosztu sum(count[j]^2) dla dopuszczonych kandydatów. Wybiera najmniejsze proxy, przy remisie najmniej miejsc L, potem wcześniejszy wariant (Natural, MD, ND). Jeśli naturalny front <=2, pomija pozostałe (brak wypełnienia i minimalny koszt strukturalny). `candidates[i].status=-1` oznacza brak próby; 0 sukces, P60_LIMIT odrzucenie budżetowe. Błędy inne niż limit przerywają analizę. `max_front` i `max_slots` dotyczą przewidywań symbolicznych; osobne limity etapu 61 kontrolują rzeczywiste fronty po pivotowaniu. AUTO nie powtarza faktoryzacji numerycznej przy niepowodzeniu. `work_bytes` to maksimum estymacji poszczególnych faz, nie suma jednoczesnych alokacji; limity nie obejmują starego faktora, wejścia, metadanych alokatora, pamięci wewnętrznej qsort ani łącznego RSS. Wybór może zwiększyć czas analizy. Nie ma gwarancji minimalnego wypełnienia ani poprawy czasu. Numeryczne łączenie superwęzłów nadal otwarte.

## Weryfikacja

`make nd`, `make nd-sanitize`, `make stage64-check`. GCC 14.2, Clang 17, ASan+UBSan. Nowe 938 testów: 800 prób na 400 losowych strukturach (ND i AUTO; 400 faktoryzacji dla pierwszych 200 przypadków), 120 układów par 2x2, 6 przypadków walidacji/limitów, 12 prób dużych układów (3 rozmiary x 4 metody, faktoryzacja dla 1600 DOF i AUTO dla 152592 DOF). Niezależna gęsta eliminacja grafowa sprawdza permutację, liczbę miejsc, front i proxy kosztu. Skalowany błąd wsteczny 7.619e-16; max błąd rozwiązania 1.554e-15. Pełna regresja 6146/6146 dla każdej konfiguracji. Log `test-results.txt` zawiera przebieg równoległego make -j4; `benchmark-serial.txt` to osobny przebieg zwykły bez równoległej kompilacji (czasy CPU pojedynczej próby, nie mediana).

## Wyniki syntetyczne (serial)

| Siatka | Metoda | slots | front | proxy | czas analizy [s] |
|---|---|---:|---:|---:|---:|
|40x40|Natural|62439|41|2602757|0.000231|
|40x40|MD|19904|63|551448|0.003397|
|40x40|ND|22329|60|620687|0.000629|
|128x128|Natural|2080895|129|269833853|0.004329|
|128x128|MD|400929|244|37100097|0.111042|
|128x128|ND|357777|192|24308715|0.008759|

AUTO wybrał MD dla 40x40, ND dla 128x128, Natural dla trójdiagonalnej 152592 DOF. Faktoryzacja 128x128 nie była wykonywana. Dla 40x40 pełna analiza+faktoryzacja: Natural 0.544790 s, MD 0.111277 s, ND 0.120651 s, AUTO 0.111674 s. Weryfikacja rozwiązania poza tym czasem. Brak rzeczywistej macierzy CCX i kodu etapów 52/54; nie wykonano benchmarku BUCKLE, MPI, ABI ani integracji z LIS.
