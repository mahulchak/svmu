#to compile the executable for svmu, the variant caller

CC = g++
CFLAGS = -g -Wall -std=c++0x 

default: svmu
svmu: svlib.o ansv.o small.o svmu.o
	$(CC) $(CFLAGS) -o svmu svlib.o ansv.o small.o svmu.o

svlib.o: svlib.cpp sv.h
	$(CC) $(CFLAGS) -c svlib.cpp

ansv.o: ansv.cpp sv.h
	$(CC) $(CFLAGS) -c ansv.cpp

small.o: small.cpp seqIO.h sv.h
	$(CC) $(CFLAGS) -c small.cpp

svmu.o: svmu.cpp sv.h seqIO.h
	$(CC) $(CFLAGS) -c svmu.cpp

clean:
	$(RM) *.o 
