g++ -O2 -ansi -pedantic -W -Wall -Wshadow -std=c++11 $(root-config --cflags) -I$PYTHIA8/include $(root-config --libs) -L$PYTHIA8/lib -lpythia8 ssbar_generate.cpp -o ssbar_generate
