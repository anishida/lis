# LIS BUCKLE – etap 69, kandydat poprawki 58

Samodzielny adapter do integracji, **nie zintegrowana biblioteka LIS/CCX**. Bazuje na źródłach prototypu etapu 68. Nie ma w dostępnym środowisku oryginalnych zmodyfikowanych źródeł LIS z etapów 52/54 ani modelu CCX 152592 DOF. Nie modyfikowano ABI LIS, nie wykonano benchmarku BUCKLE ani publikacji GitHub.

## Nowy moduł
`buckle_bridge.h`, `buckle_bridge.c`: p69_csr (real double, CSR sortowany, indeksy 0/1, pełny lub jeden trójkąt), p69_options, p69_stats, p69_factorize. Dla pełnego CSR wymagane dokładne odbicie struktury i wartości, bez automatycznego uśredniania. Duplikaty, indeksy poza zakresem, niefinitywne wartości, niezgodność wymiarów i przekroczenia limitów odrzucane. Tworzenie dolnego CSC dla K-sigma*G z sumy struktur (z zachowaniem zer strukturalnych), obliczenia przesunięcia w long double i sprawdzanie zakresu binary64. G może być NULL przy sigma=0. Zwrócone faktory rozwiązują zaokrąglony do double operator. Przy błędzie stary kontekst i statystyki pozostają niezmienione; wywołujący musi sprawdzić kod powrotu i nie użyć starego faktora jako faktora nowej macierzy. Budżet konwersji obejmuje żądane rozmiary własnych tablic, nie qsort, metadane alokatora ani inne moduły. Koszt sortowania O(m log m). Brak obsługi duplikatów, liczb zespolonych, niesymetrii, indeksów > uint32_t, rozproszonego CSR MPI. Domyślnie 1 wątek, ponieważ ostrzeżenia TSan etapu 68 są nierozstrzygnięte. OpenMP jest nadal zależnością kompilacji tej wersji.

## Schemat użycia w docelowym kodzie LIS (przykład C, bez zależności od lis.h)
```c
p56_context *ctx = NULL;
p69_options opt;
p69_defaults(&opt);
/* kp/ki/kv i gp/gi/gv: sprawdzone tablice CSR; n i nnz po sprawdzeniu zakresu. */
p69_csr k = {n, knnz, kp, ki, kv, 0, P69_FULL};
p69_csr g = {n, gnnz, gp, gi, gv, 0, P69_FULL};
int rc = p69_factorize(&ctx, &k, &g, sigma, &opt, NULL);
if (rc == 0) {
    /* scratch ma n double i nie pokrywa się z b ani x. */
    rc = p56_apply(ctx, b, x, scratch, n);
}
p56_destroy(ctx);
/* propaguj rc; nie używaj x jeśli rc != 0 */
```
W LIS trzeba zweryfikować typ i szerokość LIS_INT, indeksację CSR, zakres lokalny/globalny, typ LIS_SCALAR, konwencję znaku macierzy geometrycznej, faktyczny operator przesunięty, sposób obsługi kodów błędu oraz cykl życia preconditionera. Nie rzutować wskaźników tablic o różnych typach. Jeśli LIS już dostarcza A_sigma, przekazać ją jako K z G=NULL i sigma=0. Solver wartości własnych pozostaje po stronie LIS/CCX; adapter nie liczy modów wyboczenia. Dla MPI konieczna osobna integracja i walidacja; nie wolno interpretować lokalnego bloku jako całej macierzy. Nie aktywować OpenMP produkcyjnie przed weryfikacją TSan. Klucz cache musi uwzględniać operator i przesunięcie oraz konfigurację; adapter nie włącza automatycznego cache.

## Testy odtwarzalne
```
make bridge
make bridge-sanitize
make stage69-check
```
GCC 14.2; nowe 436/436 przypadków przechodzi w obu konfiguracjach. 360 losowych problemów z różnymi formatami, bazami i metodami porządkowania, 60 problemów z wymuszonymi pivotami blokowymi po przesunięciu, 15 przypadków walidacji (w tym poprawna budowa początkowa), 1 test dużego układu. Błąd max x 5.551e-16, komponentowy błąd wsteczny 4.051e-16. Pierwsza kompilacja sanitizera wykryła ostrzeżenie GCC o możliwie niezainicjalizowanym buforze testowym, poprawiono przez calloc, obie końcowe kompilacje i testy PASS. Log bridge-results.txt zawiera końcowe wyniki. Ponowna regresja stage68-check przerwana przez limit czasu środowiska, zakończone testy sanitizowane etapów 55–64 (6146 przypadków) PASS, nie zaliczać pozostałych jako PASS na podstawie tego przebiegu. Kod solvera bez zmian względem etapu 68. Wcześniejszy wynik 8857/8857 dotyczy etapu 68, nie nowego pełnego przebiegu. Log regression-stage68.txt zawiera przerwany przebieg. Nie zweryfikowano TSan na innym runtime.

## Pomiar syntetyczny
152592 DOF, K pełne trójdiagonalne 457774 wpisów, G diagonalna 152592 wpisy; sigma=0.5, 305183 wpisy dolnego CSC. Czas ścienny konwersji, porządkowania, analizy i faktoryzacji 0.212835 s, pojedyncze zastosowanie 0.002509 s. Bufory konwersji 17700644 B, monitorowana pamięć numeryczna 23350680 B; wartości z różnych faz, nie dodawać jako zmierzonego RSS. Pojedynczy pomiar bez sanitizerów, nie benchmark CCX ani gwarancja wydajności. ldd: libm, libgomp, libc, loader systemowy, bez BLAS/LAPACK.

## Warunki etapu 70
Wymagane: źródła LIS z poprawkami 52/54 i źródła/konfiguracja CCX, rzeczywisty plik .inp oraz wszystkie INCLUDE (lub eksport dokładnych macierzy, RHS i konfiguracji). Sprawdzić znak i przesunięcie, porównać z etapem 52, przebudować LIS/CCX, wykonać pełną regresję, pomiary wall time i peak RSS, certyfikację sześciu modów z residuum oryginalnego problemu i testy MPI/ABI. Bez tych danych nie ma możliwości ukończenia wydania produkcyjnego.
