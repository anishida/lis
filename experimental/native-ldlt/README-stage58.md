# LIS BUCKLE — etap 58 / kandydat poprawki 47

Samodzielny prototyp C11, oparty na rzeczywistym archiwum etapu 57. Nie jest poprawką do oryginalnego LIS; źródła etapów 52/54 i macierz CCX nie były dostępne. Zachowano format cache i API etapów 55–57.

## Algorytm
`p58_factorize` przyjmuje pełną, skończoną, dokładnie symetryczną macierz double w układzie row-major, symboliczne panele o kolejnych zakresach kolumn i posortowane listy wierszy zewnętrznych. Kolejność wejścia jest początkową kolejnością DOF. `max_n` ogranicza rozmiar gęstej referencji (obowiązkowy). `tau` należy do (0,1].

1. Akceptuj niezerowy pivot skalarny, jeśli |a_kk| >= tau * max_{i>k}|a_ik|.
2. W przeciwnym przypadku szukaj niezerowego partnera j w tym samym panelu; wybierz największe |a_jk| spośród dopuszczalnych kandydatów. Sprawdź skalowany wyznacznik (próg 8 LDBL_EPSILON względem sumy modułów składników) i każdy mnożnik bloku (moduł <= 1/tau i reprezentowalny jako double).
3. Symetrycznie zamień j z k+1 w aktywnym Schurze, zamień odpowiednie wiersze wcześniej obliczonego L i permutację. Zamiana jest dopuszczalna tylko przy jednakowej przynależności obu wierszy do list zewnętrznych wszystkich wcześniejszych paneli.
4. Aktualizuj symetryczny Schur z użyciem long double / fmal; zapisuj do double po kontroli zakresu. W L miejsce pod przekątną wewnątrz bloku D jest zerem.
5. Sprawdź, czy poza zadanym wzorcem L nie ma żadnych niezerowych wartości (również drobnych pozostałości zaokrągleń). Zapakuj przez p56_rebuild, transakcyjnie zastępując kontekst. Stats zmieniane wyłącznie przy sukcesie.

`P58_NO_PIVOT` oznacza brak dopuszczalnego pivota, a nie dowód osobliwości. `P58_PATTERN` oznacza niezgodność wzorca symbolicznego. Pozostałe kody pochodzą z P55. Brak automatycznego wyboru alternatywnego pivota 1x1, brak pełnego Bunch–Kaufmana i brak opóźnionego pivotowania. Kryterium lokalne nie jest dowodem globalnej stabilności. Wzorzec nie jest automatycznie obliczany. Niezerowe elementy poza wzorcem są odrzucane, nie odcinane. Zakresy i precyzja long double zależą od platformy.

## Testy
`make stage58-check` uruchamia wszystkie testy optymalizowane i ASan+UBSan (CC i SAN_CC można nadpisać). `make clean` usuwa binaria.
- 314 testów pakowania, 318 adaptera, 560 cache (niezmienione).
- 300 generowanych układów 10x10 (2 skalarne i 4 blokowe pivoty), 30 zapisów/odczytów cache, 120 układów 8x8 z niezerowymi aktualizacjami Schura wewnątrz oraz między panelami i 11 testów szczególnych.
- Niezależna rekonstrukcja P A P^T z rozpakowanych L,D oraz kontrola rozwiązania względem oryginalnej macierzy. Skale 1e-100, 1, 1e100.
- Wyniki: 461/461 nowych; 1653/1653 łącznie w każdej konfiguracji. W 420 generowanych układach 1680 par i 955 zamian (bez dodatkowych testów szczególnych). Maks. błąd bezwzględny 4.441e-16, skalowany błąd wsteczny 1.520e-16 w nowych testach. Definicja błędu wstecznego: ||Ax-b||_inf/(||A||_inf||x||_inf+||b||_inf).
- `test-results.txt` zawiera rzeczywiste logi kompilacji, uruchomienia i ldd. Jedynie libc/libm i loader, bez BLAS/LAPACK.

## Ograniczenia i dalsza integracja
To gęsta referencja O(n^3) dla eliminacji i O(n^2) pamięci (same dwie tablice double: 16*n*n bajtów; dodatkowo wejście i bufor pakowania, m.in.). Wyszukiwanie partnerów może zwiększać koszt. Nie używać dla 152592 DOF; samych dwóch gęstych tablic potrzeba ok. 347 GiB. Docelowo przenieść logikę do rzadkich paneli i ich aktualizacji. Format cache zachowano; klucz musi rozróżniać algorytm, tau, macierz, wzorzec i kolejność. Nie wykonano CCX, MPI, pełnego make check LIS ani weryfikacji ABI. Nie publikowano na GitHubie. Etap 52 nadal jest zalecany dla macierzy CCX na podstawie przekazanego raportu, nie nowych pomiarów.
