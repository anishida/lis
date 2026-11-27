# LIS BUCKLE — etap 57 / kandydat poprawki 46

Samodzielny prototyp C11/POSIX oparty na etapie 56. Brak BLAS, OpenBLAS i LAPACK; libm i libc. Nie jest zintegrowaną poprawką LIS: nadal brak źródeł etapów 52/54 i macierzy CCX. Prototyp przyjmuje gotowe faktory i nie realizuje eliminacji ani wyboru pivotów.

## API

`p57_save(ctx, path, key)` zapisuje niezmienny kontekst `p56_context` do pliku tymczasowego w tym samym katalogu (mkstemp, uprawnienia 0600), wykonuje fflush, fsync pliku, fclose i rename. Nieudany zapis usuwa plik tymczasowy; nie zastępuje poprzedniego docelowego pliku przed sukcesem. Nie ma fsync katalogu, więc nie gwarantuje trwałości wpisu katalogowego po utracie zasilania. Nie ma blokady między piszącymi; ostatni udany rename wygrywa. Czytanie i zapis wymagają kontrolowanego katalogu i synchronizacji życia kontekstu po stronie wywołującego.

`p57_load(&ctx, path, key, max_bytes, max_n)` czyta regularny plik, ogranicza jego rozmiar oraz liczbę niewiadomych, sprawdza format i klucz, waliduje liczniki i rozmiary przed alokacjami, dekoduje dane, sprawdza CRC, następnie wywołuje istniejącą walidację `p56_rebuild` / `p55_pack`. Dopiero po sukcesie i zamknięciu pliku zastępuje kontekst. W razie błędu stary kontekst niezmieniony. Limity pliku i n nie są ścisłym limitem RSS: dekodowanie przejściowo używa nieskompresowanych indeksów i kopii faktorów. `P57_MISS` oznacza brak pliku albo niezgodny identyfikator/klucz; `P57_IO` błąd I/O; `P57_CORRUPT` złą strukturę lub checksum; `P57_LIMIT` przekroczone limity. Istniejące kody P55_* zachowano.

Klucz 32-bajtowy dostarcza wywołujący. Prototyp nie oblicza klucza macierzy. Integracja musi powiązać go np. z SHA-256 kanonicznej struktury i wartości macierzy, przesunięcia spektralnego, permutacji, ustawień, wersji algorytmu/kompilacji i rozkładu MPI. CRC32 wykrywa przypadkowe uszkodzenia; nie jest uwierzytelnieniem ani ochroną przed celową podmianą. Nie należy ładować cache z niezaufanego źródła. Nie sprawdzono interoperacyjności na innych architekturach. `long double` zależy od platformy; po odczycie normalizacja D jest przeliczana. Dodatkowe pole original_b zachowuje dokładny współczynnik bloku.

## Format LDLTC001
Wszystkie liczby całkowite bez znaku w little-endian, double jako 8 bajtów binary64 little-endian. Kolejno:
- 8 bajtów magic, 32 bajty key, cztery uint32: n, np, nb, flaga permutacji (0/1);
- n double diagonali;
- nb par (uint32 k, double off[k]), uporządkowane nieprzecinające się bloki;
- opcjonalne n uint32 perm[k] (factor-to-original);
- np paneli: uint32 width, nrows, index_bytes_per_entry (2 lub 4); nrows indeksów (uint16 względne do końca panelu albo uint32 absolutne); width*(width-1)/2 + width*nrows double w kolejności kolumnowej etapu 55. first wynika z sumy width poprzednich paneli;
- uint32 CRC-32/ISO-HDLC całego poprzedniego strumienia, bez trailing bytes.

Zapis nie zawiera wskaźników ani reprezentacji long double. Rozmiar pliku: 60 + 8n + 12nb + 4n*flag + suma(12 + nrows*index_bytes_per_entry + 8*nvalues) bajtów. Minimalny przykładowy test 2x2 z permutacją ma 116 B. Czas I/O, hit/miss oraz RSS nie zostały zbenchmarkowane. Brak gwarancji poprawności cache dla innej macierzy, jeśli integrator błędnie poda ten sam klucz.

## Testy
`make full-check` — 314 starych testów pakowania + 318 adaptera + 560 nowych cache, osobno w O3 i ASan+UBSan. 300 nowych przypadków opiera się na deterministycznych losowych faktorach, niezależnej rekonstrukcji A i porównaniu rozwiązania z przed/po cache za pomocą memcmp; 116 modyfikacji najmłodszego bitu po jednej na bajt pliku, 116 obcięć, 15 celowo błędnych struktur z ponownie wyliczonym CRC oraz pozostałe przypadki, w tym limity, brak pliku, FIFO, rollback i 65538 DOF z indeksami 16/32-bit. Wszystkie zakończyły się sukcesem w lokalnym środowisku; logi w test-results.txt. Nie testowano utraty zasilania, awarii alokacji, konkurencyjnych zapisów, MPI, innych platform ani pełnej biblioteki LIS.

## Kolejny etap
Dynamiczna symetryczna zamiana pivotów podczas eliminacji w obrębie ścisłego superwęzła. Potem kryteria Bunch–Kaufmana i opóźnione pivotowanie. Integracja z LIS nadal wymaga źródeł i danych CCX; ABI publiczne LIS nie było badane. Nie publikowano zmian na GitHubie.
