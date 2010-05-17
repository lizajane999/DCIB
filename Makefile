# Make file for dcib
# Lisa Miller 2010

dcib: dcib.c histo2.c utils.c utils_ps.c utils_IB.c
	gcc -g -Wall -O -lm -o dcib dcib.c histo2.c utils.c utils_ps.c utils_IB.c -I.
	
tidy:
	rm *.o

clean:
	rm -f dcib
