g++ -std=c++17 -O2 -fPIC -shared -I ./glm -o main.so main.cpp
g++ -std=c++17 -O2 -fPIC -shared -I ./glm -o spectral.so spectral2.cpp
g++ -std=c++17 -O2 -fPIC -shared -I ./glm -o lindstrom_turk.so lindstrom_turk.cpp