# LIS BUCKLE – etap 60 / kandydat poprawki 49

Samodzielny prototyp C11 rozwijający etap 59. Nie jest poprawką zastosowaną do rzeczywistych źródeł LIS. Bez BLAS, OpenBLAS i LAPACK.

## Algorytm
`p60_factorize` używa klasycznego dolnego Bunch–Kaufmana z etapu 59 na gęstej macierzy Schura. Oryginalne panele symboliczne są uporządkowane liniowo. Gdy wybrany pivot wymaga niedostępnej kolumny, bieżąca kolumna jest zamieniana z ostatnim dostępnym kandydatem i trafia do sufiksu oczekującego. Po wyczerpaniu kandydatów kolejny panel obejmuje sufiks i nowe kolumny. Zamiany obejmują poprzednie kolumny L i permutację. W ostatnim panelu wszystkie pozostałe kolumny są dostępne. Zero w całej aktywnej kolumnie oznacza odrzucenie.

Po eliminacji odtwarzane są panele i ich listy wierszy na podstawie faktycznych niezerowych współczynników L w końcowej kolejności faktorów; nie usuwamy małych niezerowych elementów. Panele bez wyeliminowanych kolumn nie są zapisywane. Blok 2x2 nie jest dzielony przez granicę wyjściowego panelu. Format pozostaje zgodny z adapterem i cache z etapów 56–57. Wejściowy wzorzec jest sprawdzany względem niezerowych elementów A. To nie jest symboliczna analiza rzadkiej eliminacji.

## API i limity
Patrz `delayed_ldlt.h`. `max_n` ogranicza rozmiar gęstej macierzy, `max_front` szerokość obszaru kandydatów wraz z opóźnionymi kolumnami na wejściu do panelu, `max_delays` całkowitą liczbę zdarzeń opóźnienia (kolumna może być opóźniona wielokrotnie), `max_pending` maksymalną liczbę kolumn w sufiksie oczekującym, a `max_extra_slots` wzrost liczby przechowywanych double w L względem pierwotnego wzorca paneli (uwzględnia zera w pakowanych panelach). To nie jest liczba nowych pozycji strukturalnych ani całkowity limit RAM. Budżet miejsc sprawdzany po eliminacji i przed alokacją pakowanych wartości. Limit `max_growth=0` wyłącza kontrolę wzrostu, inna wartość musi być >=1. Przekroczenie limitu frontu, opóźnień, sufiksu lub miejsc: `P60_LIMIT`. Pozostałe statusy zgodne z etapami wcześniejszymi. Kontekst i statystyki zmieniają się tylko przy sukcesie. Wywołujący serializuje przebudowę względem apply.

## Testy
`make stage60-check` – wszystkie zestawy zwykłe i ASan/UBSan. `make delayed` i `make delayed-sanitize` – wyłącznie nowe testy. `make clean` – usuwa programy testowe. Wykonano GCC 14.2 (-O3) i Clang 17 (-O1 z ASan/UBSan). Patrz `test-results.txt`.

- 2865/2865 asercji testowych w każdej konfiguracji, w tym 564/564 nowych (400 losowych macierzy wielopanelowych, 120 rzadkich macierzy z odległymi parami, 31 roundtrip cache, 13 przypadków szczególnych).
- W 400 losowych przypadkach etap 59 rozwiązał 40, etap 60 wszystkie 400; zestaw celowo obciąża ograniczenie panelowe. Skale 1e-100, 1, 1e100. Rekonstrukcja PAP^T i residuum badane niezależnie.
- W 520 generowanych przypadkach 6236 opóźnień, 1227 bloków, 5196 zamian, maksymalny front 10, maksymalny sufiks 9. Maksymalny wzrost 13.00686. Maksymalny błąd bezwzględny rozwiązania 2.782e-13; maksymalny skalowany błąd wsteczny 3.238e-16 (także przypadki szczególne i cache).
- Limit dodatkowych miejsc: przykład 3x3 ma 2 pierwotne miejsca i 3 miejsca po pakowaniu. Budżet 0 odrzuca, budżet 1 akceptuje. Potwierdzono zachowanie wcześniejszego faktora po błędzie.
- `ldd` programu test_delayed: libm, libc i loader systemowy; nie ma BLAS/LAPACK.

## Ograniczenia i kolejny krok
To model sekwencyjnych paneli, nie algorytm multifrontalny ani drzewo eliminacji. Dwie gęste macierze robocze wymagają co najmniej 16*n*n bajtów (bez A wejściowego, faktorów i buforów pomocniczych). Nie ma gwarancji znalezienia pivota przy skończonych limitach. Nie ma pełnego limitu bajtów ani pomiaru RSS. Wzrost miejsc L nie jest kontrolą liczby nowo powstałych pozycji względem konkretnego symbolicznego wzorca. Dla CCX 152592 DOF podejście gęste nie jest praktyczne. Brak źródeł etapów 52/54 oraz rzeczywistej macierzy CCX uniemożliwia integrację, benchmark, MPI i ABI. Nie publikowano na GitHubie. Etap 61: przenieść mechanizm do rzadkich frontów z drzewem eliminacji, rzeczywistym budżetem RAM i aktualizacjami blokowymi; dopiero później SIMD.
