BIN = ../../bin
LOGDIR=../../log
LIBS=-lcrlibm

SEQA=GCGCGTGCGCGGAAGGAGCCAAGGTGAAGTTGTAGCAGTGTGTCAGAAGAGGTGCGTGGCACCATGCTGTCCCCCGAGGCGGAGCGGGTGCTGCGGTACCTGGTCGAAGTAGAGGAGTTG
SEQB=GACTTGTGGAACCTACTTCCTGAAAATAACCTTCTGTCCTCCGAGCTCTCCGCACCCGTGGATGACCTGCTCCCGTACACAGATGTTGCCACCTGGCTGGATGAATGTCCGAATGAAGCG

CC=g++
DEFS=
GLOBAL_DEFS=
WARNINGS=-Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wmissing-declarations -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized
CFLAGS+= $(WARNINGS) -O3 $(DEFS) $(GLOBAL_DEFS) -std=c++11 -march=native

TARGET=$(BIN)/TKFLOG_CACHING_ROUND_UP

DEPS=TKF91Sequential.h
OBJS=TKF91Sequential.o Main.o
OUT_OBJS=$(addprefix $(BIN)/,$(OBJS))

valgrind: CFLAGS+= -g
clang: CC=clang++
clang: WARNINGS=-Weverything

all:  $(TARGET)

clang:  $(TARGET)

cmd1: $(TARGET)
	$(TARGET) --sequence-1 ACATA --sequence-2 CAATT --lambda 1 --mu 2 --tau 0.1 --pa 0.27 --pc 0.24 --pg 0.26 --pt 0.23

cmd3: $(TARGET)
	$(TARGET) --sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA \
	--sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA --lambda 1 --mu 2 \
	--tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25 -v 

callgrind: $(TARGET)
	mkdir -p $(LOGDIR)
	valgrind --tool=callgrind --callgrind-out-file=$(LOGDIR)/callgrind.output.%p \
	$(TARGET) --sequence-1 $(SEQA) --sequence-2 $(SEQB) --lambda 1 --mu 2 \
	--tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25

bench: $(TARGET)
	cd .. && ../bin/TKFLOG_CACHING_ROUND_UP -t

valgrind: $(TARGET)
	valgrind $(TARGET) --sequence-1 $(SEQA) --sequence-2 $(SEQB) --lambda 1 --mu 2 \
	--tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.2

vgdb: $(TARGET)
	valgrind --vgdb=yes --vgdb-error=0 $(TARGET) --sequence-1 $(SEQA) --sequence-2 $(SEQB) --lambda 1 --mu 2 \
	--tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25

$(TARGET): $(OUT_OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(OUT_OBJS): $(BIN)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean all profile 

clean:
	rm -f $(BIN)/*.o
	rm -f $(TARGET)
