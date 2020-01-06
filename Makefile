CXXFLAGS+=-std=c++11 -O3 -fno-tree-vectorize -DNDEBUG -Wall -Wextra -g
TARGET=memstress

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	$(RM) $(TARGET) *.o
