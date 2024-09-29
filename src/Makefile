CC=g++
CFLAGS= -std=c++17 -O3 -funroll-loops -fomit-frame-pointer 

#TARGET = star-cdp
#SRC = $(#wildcard *.cpp)
#OBJ = $(#patsubst %.cpp, %.o, $(#SRC))
all: star-cdp

#$(#TARGET): $(#OBJ)
#	$(#CC) -o $@ $^

#%.o: %.cpp
#	$(#CC) $(#CFLAGS) -I boost_1_80_0 $< -o $@

star-cdp: driver.cpp 
	$(CC) -DUSE_SHRT $(CFLAGS) -I boost_1_80_0 -o star-cdp driver.cpp 


.PHONY: clean
clean:
	rm star-cdp
