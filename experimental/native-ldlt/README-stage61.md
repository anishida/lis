# LIS BUCKLE – etap 61 / kandydat poprawki 50

Samodzielny prototyp C11. Kontynuacja etapów 55–60. Nie jest poprawką
zintegrowaną z oryginalnym LIS. Bez BLAS, LAPACK i OpenBLAS.

## Algorytm

`p61_factorize` przyjmuje dolny trójkąt symetrycznej macierzy w CSC,
indeksy 0-based, wiersze posortowane rosnąco bez duplikatów, brakująca
diagonala oznacza zero. `col` ma n+1 elementów, `row` i `val` po col[n]
elementów. Wywołujący dostarcza tablicę `parent` o n elementach: rodzic
ma większy indeks, korzeń to UINT32_MAX. Każdy strukturalny wpis
pozadiagonalny musi łączyć kolumnę z jej przodkiem. Sprawdzamy warunek
za pomocą iteracyjnego DFS i przedziałów wejścia/wyjścia. Las jest
obsługiwany. Kolejność rodzic > dziecko wystarcza do obliczeń bottom-up.

Dla każdego frontu: zbieramy kolumnę własną, opóźnione kolumny dzieci,
indeksy oryginalnych wpisów CSC i indeksy aktualizacji dzieci. Montujemy
lokalną gęstą macierz, sumując wkłady Schura. Wybieramy pivot 1×1 lub
2×2 ze strategią Bunch–Kaufmana ograniczoną do aktywnych, w pełni
zmontowanych kandydatów. Gdy wybrany partner nie jest aktywny,
przenosimy bieżącego kandydata do aktualizacji rodzica. W korzeniu
brak dopuszczalnego pivota powoduje błąd. Wiersze aktualizacyjne nie
muszą być w pełni zmontowane: ich częściowe maksima nie są globalnymi
maksimami; strategia nie jest tożsama z globalnym Bunch–Kaufmanem.
Aktualizacje dzieci są zwalniane po zmontowaniu frontu nadrzędnego.

Mnożniki zapisujemy wg oryginalnych indeksów, po zakończeniu eliminacji
mapujemy do finalnej permutacji i sortujemy. Każdy pivot tworzy panel
o szerokości 1 lub 2; nie ma jeszcze łączenia ścisłych superwęzłów.
Używamy dotychczasowego p56_rebuild, p56_apply i p57_save/load.
Poprzedni kontekst i statystyki pozostają niezmienione przy błędzie.

## Limity

- max_front: całkowita lokalna szerokość, wliczając wiersze aktualizacji.
- max_work_bytes: sumaryczny limit żądanych bajtów alokacji TEGO modułu,
  obejmuje fronty, aktualizacje, metadane i tymczasowe pakowanie.
  Nie obejmuje metadanych malloc, stosu/qsort, wejścia CSC, starego
  kontekstu ani alokacji końcowego p56_rebuild. Nie jest limitem RSS.
- max_factor_slots: limit liczby wartości double w pakowanym L;
  obejmuje też strukturalne zero wewnątrz pivota 2×2.
- max_delays: liczba przeniesień (ta sama kolumna może być liczona wielokrotnie).
- max_growth: 0 = brak limitu; w przeciwnym razie >=1. Wzrost odnosi się
  do obserwowanych wartości lokalnych frontów, w tym częściowych sum
  montażu, a nie globalnego wzrostu klasycznego Bunch–Kaufmana.

Limit pamięci jest sprawdzany przed każdą alokacją własną. API zakłada
poprawne rozmiary rzeczywiście dostępnych buforów (jak inne API C).
Serializacja równoczesnej przebudowy i odczytu pozostaje obowiązkiem
wywołującego. Obliczenia są sekwencyjne. Nie ma globalnej macierzy n*n,
ale lokalny front może być duży; ma koszt pamięci O(m^2) i aktualizacje
O(m^2) na pivot. Rozmiary wkładów i faktorów zależą od fill-in i drzewa.

## Kompilacja i regresja

```
make stage61-check
make multifrontal
make multifrontal-sanitize
```

543 nowych testów: 300 losowych drzew (macierze diagonalnie dominujące,
symetryczne nieokreślone), 200 układów z odległymi blokami, 3 dodatkowe
przypadki strukturalne, 10 testów błędów/cache, 26 limitów pamięci,
4 większe układy. 2865 poprzednich testów ponownie uruchomiono.
Razem 3408 zwykłych i 3408 z ASan/UBSan. Zobacz `test-results.txt`.
Nowy zestaw: max błąd rozwiązania 8.549e-15, max komponentowy skalowany
błąd wsteczny 4.831e-16. Skale 1e-100..1e100. 78657 bloków 2x2,
78659 przeniesień (liczniki udanych przypadków run()).

Duże syntetyczne przypadki: trójdiagonalny 152592 DOF, 305183 wpisów
CSC, front 2, 152591 miejsc w L, peak_work_bytes=22737252;
nieokreślony blokowo-diagonalny 152592 DOF, 76296 pivotów 2x2,
76296 przeniesień, peak_work_bytes=20906128. Siatka 32x32: 1024 DOF,
front 33, 31775 miejsc w L. Nie są to macierze CCX. Czasy CPU w logu
nie są czasami pełnego BUCKLE; zawierają odpowiednio faktoryzację z
pakowaniem i oddzielnie rozwiązanie wraz z weryfikacją.

## Następny krok

Automatyczna analiza symboliczna CSC, rzeczywiste drzewo eliminacji,
łączenie superwęzłów, optymalizacja planowania drzewa, testy większych
macierzy nieregularnych oraz budżet obejmujący końcowe kopie faktorów.
Integracja z oryginalnym LIS wymaga źródeł etapów 52 i 54, a benchmark
CCX wymaga macierzy. Nie wykonano testów całego LIS, MPI ani ABI.
