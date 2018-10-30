CC=g++
SHARKLIB=/usr/lib/
SHARKINC=/usr/include/
LDLIBS=-lshark
LDFLAGS=-L${SHARKLIB} -Wl,-rpath,${SHARKLIB}
#CXXFLAGS=-Wall -pedantic -I${SHARKINC} -Ofast
CXXFLAGS=-Wall -pedantic -I${SHARKINC} -O4
#CXXFLAGS=-Wall -pedantic -I${SHARKINC} -ggdb -DDEBUG

OBJECTS=main.o RunParameter.o Benchmarks.o CCVIL.o Archive.o\
				F1.o F2.o F3.o F4.o F5.o F6.o F7.o F8.o F9.o F10.o\
				F11.o F12.o F13.o F14.o F15.o F16.o F17.o F18.o F19.o F20.o 

CCVIL: $(OBJECTS)
	$(CC) $(CXXFLAGS) -o CCVIL $(OBJECTS) $(LDFLAGS) $(LDLIBS)

main.o: main.cpp Header.h benchmark/RunParameter.h benchmark/Benchmarks.h CCVIL.h Archive.h\
	benchmark/F1.h benchmark/F2.h benchmark/F3.h benchmark/F4.h benchmark/F5.h benchmark/F6.h benchmark/F7.h benchmark/F8.h benchmark/F9.h benchmark/F10.h\
	benchmark/F11.h benchmark/F12.h benchmark/F13.h benchmark/F14.h benchmark/F15.h benchmark/F16.h benchmark/F17.h benchmark/F18.h benchmark/F19.h benchmark/F20.h
	$(CC) $(CXXFLAGS) -c main.cpp 

CCVIL.o: CCVIL.h benchmark/Benchmarks.h Archive.h CCVIL.cpp
	$(CC) $(CXXFLAGS) -c CCVIL.cpp

Benchmarks.o: benchmark/RunParameter.h benchmark/Benchmarks.h benchmark/Benchmarks.cpp
	$(CC) $(CXXFLAGS) -c benchmark/Benchmarks.cpp

RunParameter.o: benchmark/RunParameter.h benchmark/RunParameter.cpp
	$(CC) $(CXXFLAGS) -c benchmark/RunParameter.cpp

Archive.o: Archive.h Archive.cpp
	$(CC) $(CXXFLAGS) -c Archive.cpp

F1.o: benchmark/F1.h benchmark/Benchmarks.h benchmark/F1.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F1.cpp

F2.o: benchmark/F2.h benchmark/Benchmarks.h benchmark/F2.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F2.cpp

F3.o: benchmark/F3.h benchmark/Benchmarks.h benchmark/F3.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F3.cpp

F4.o: benchmark/F4.h benchmark/Benchmarks.h benchmark/F4.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F4.cpp

F5.o: benchmark/F5.h benchmark/Benchmarks.h benchmark/F5.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F5.cpp

F6.o: benchmark/F6.h benchmark/Benchmarks.h benchmark/F6.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F6.cpp

F7.o: benchmark/F7.h benchmark/Benchmarks.h benchmark/F7.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F7.cpp

F8.o: benchmark/F8.h benchmark/Benchmarks.h benchmark/F8.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F8.cpp

F9.o: benchmark/F9.h benchmark/Benchmarks.h benchmark/F9.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F9.cpp

F10.o: benchmark/F10.h benchmark/Benchmarks.h benchmark/F10.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F10.cpp

F11.o: benchmark/F11.h benchmark/Benchmarks.h benchmark/F11.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F11.cpp

F12.o: benchmark/F12.h benchmark/Benchmarks.h benchmark/F12.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F12.cpp

F13.o: benchmark/F13.h benchmark/Benchmarks.h benchmark/F13.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F13.cpp

F14.o: benchmark/F14.h benchmark/Benchmarks.h benchmark/F14.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F14.cpp

F15.o: benchmark/F15.h benchmark/Benchmarks.h benchmark/F15.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F15.cpp

F16.o: benchmark/F16.h benchmark/Benchmarks.h benchmark/F16.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F16.cpp

F17.o: benchmark/F17.h benchmark/Benchmarks.h benchmark/F17.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F17.cpp

F18.o: benchmark/F18.h benchmark/Benchmarks.h benchmark/F18.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F18.cpp

F19.o: benchmark/F19.h benchmark/Benchmarks.h benchmark/F19.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F19.cpp

F20.o: benchmark/F20.h benchmark/Benchmarks.h benchmark/F20.cpp
	$(CC) $(CXXFLAGS) -c benchmark/F20.cpp

.PHONY : clean clrout
clean:
	rm -f CCVIL $(OBJECTS)

clrout:
	rm -fr result trace out outout
