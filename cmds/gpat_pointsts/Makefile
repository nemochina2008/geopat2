
PROG = gpat_pointsts

LIBS = \
../../lib/argtable/libargtable3.a \
../../lib/tools/libtools.a

ifndef PREFIX
	PREFIX = /usr/local
endif

all: $(PROG)

$(PROG): main.o 
	$(CC) -o $(PROG) main.o $(LIBS) $(CFLAGS) $(EXTFLAGS)

main.o: main.c ../../lib/argtable/argtable3.h
	$(CC) $(CFLAGS) $(EXTFLAGS) -c main.c

install: $(PROG)
	mv -f $(PROG) $(PREFIX)/bin

clean:
	rm -f *.o $(PROG)