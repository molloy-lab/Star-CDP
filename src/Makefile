CC = g++
CFLAGS = -std=c++17 -O3 -funroll-loops -fomit-frame-pointer


UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    BOOST_INCLUDE_PATH = /usr/include
    BOOST_LIB_PATH = /usr/lib
else ifeq ($(UNAME_S),Darwin)  # macOS
    BOOST_INCLUDE_PATH = /opt/homebrew/include
    BOOST_LIB_PATH = /opt/homebrew/lib
endif


TARGET = star-cdp
SRC = driver.cpp


$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -I $(BOOST_INCLUDE_PATH) -L $(BOOST_LIB_PATH) -lboost_system -o $(TARGET) $(SRC)

.PHONY: clean
clean:
	rm -f $(TARGET)
