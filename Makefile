#
# 2010/12/16 by D. Grant Starkweather
#

MVPVERSION = 0.0.0

HFLS	= mvptree.h
OBJS	= mvptree.o

CC	= cc

CFLAGS	= -g -O3 -I. -Wall $(DEFINES)

CPPFLAGS += -pthread -I /usr/local/include

LDFLAGS	=
RANLIB	= ranlib

DESTDIR	= /usr/local
TEST	= testmvp
TEST2	= testmvp2
TEST3   = imget

LIBRARY	= libmvptree.a

DEPS_LIBS = -lm

IMGET_DEFINES = -Dcimg_use_xshm -Dcimg_use_xrandr -Dcimg_use_jpeg
IMGET_LIBS = -lX11 -lXext -lXrandr -ljpeg

#uncomment to use libprng.a library for random number generation
#DEPS_LIBS += /usr/local/lib/libprng.a


all : $(TEST) $(TEST2)

clean :
	rm -f a.out core *.o *.t
	rm -f $(LIBRARY) $(UTIL) $(TEST) $(TEST2) $(TEST3)

install : $(HFLS) $(LIBRARY) 
	install -c -m 444 $(HFLS) $(DESTDIR)/include
	install -c -m 444 $(LIBRARY) $(DESTDIR)/lib
	$(RANLIB) $(DESTDIR)/lib/$(LIBRARY)

$(LIBRARY) : $(OBJS) $(HFLS)
	ar cr $(LIBRARY) $?
	$(RANLIB) $@

tests : $(TEST) $(TEST2) $(TEST3)

$(TEST) : $(LIBRARY) $(TEST).o 
	rm -f $@
	$(CC) $(CFLAGS) $(LDFLAGS) $(TEST).o $(LIBRARY) $(DEPS_LIBS)
	mv a.out $@

$(TEST2): $(LIBRARY) $(TEST2).o
	rm -f $@
	$(CC) $(CFLAGS) $(LDFLAGS) $(TEST2).o $(LIBRARY) $(DEPS_LIBS)
	mv a.out $@

imget: $(LIBRARY) imget.o
	rm -f $@
	g++ $(CFLAGS) $(IMGET_DEFINES) $(CPPFLAGS) $(LDFLAGS) imget.o $(LIBRARY) $(DEPS_LIBS) $(IMGET_LIBS)
	mv a.out $@

.c.o :
	rm -f $@
	$(CC) $(CFLAGS) -c $< -o $@

imget.o : imget.cpp
	rm -f $@
	g++ $(CFLAGS) $(CPPFLAGS) -c imget.cpp -o $@
