all: density_split

density_split:
	make -C src -f Makefile install

clean:
	make -C src -f Makefile clean
	rm -f bin/*.exe
