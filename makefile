SRCS = $(shell find . -maxdepth 1 -name "*.c*") 

OBJS = $(addsuffix .o, $(basename $(SRCS)))

EXEC = hybrid_cec

LIBS = -g -lkissat -LhKis/build
LIBS = -lkissat -LhKis/build

CXXFLAGS=-Ofast -IhKis -mavx2 -std=c++17 -march=native

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC)
