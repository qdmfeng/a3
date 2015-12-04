CC = g++
CFLAGS = -g -O2 
LIBS = -lm

all: chamber

chamber:  advanced_main.o raytracer.o util.o light_source.o scene_object.o bmp_io.o
	$(CC) $(CFLAGS) -o chamber \
	advanced_main.o raytracer.o util.o light_source.o scene_object.o bmp_io.o $(LIBS)

raytracer:  main.o raytracer.o util.o light_source.o scene_object.o bmp_io.o
	$(CC) $(CFLAGS) -o raytracer \
	main.o raytracer.o util.o light_source.o scene_object.o bmp_io.o $(LIBS)

clean:
	-rm -f core *.o
	-rm raytracer
	-rm chamber
	



