all: smudge_phenos

smudge_phenos: smudge_phenos.cc
	g++ -Ofast -funroll-loops -I. -Wall -std=c++0x smudge_phenos.cc -o smudge_phenos -lboost_thread -lboost_system -lpthread
