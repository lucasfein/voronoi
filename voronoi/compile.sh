g++ --std=c++17 -o "${BASH_SOURCE%/*}/voronoi.exe" "${BASH_SOURCE%/*}/voronoi.cc" -O3 -funroll-loops -lboost_filesystem -lboost_system -Wall -Wextra -Werror -pedantic