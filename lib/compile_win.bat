@echo off
g++ -std=c++17 -O2 -shared -I./glm -o main.dll main.cpp -static-libgcc -static-libstdc++
g++ -std=c++17 -O2 -shared -I./glm -o lindstrom_turk.dll lindstrom_turk.cpp -static-libgcc -static-libstdc++
g++ -std=c++17 -O2 -shared -I./glm -o spectral.dll spectral.cpp -static-libgcc -static-libstdc++
g++ -std=c++17 -O2 -shared -I./glm -o structure_aware.dll structure_aware.cpp -static-libgcc -static-libstdc++
pause
