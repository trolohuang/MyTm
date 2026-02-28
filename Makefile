.PHONY: clean

GSL_DIR ?=
GSL_LIB ?=
MKL_DIR ?=
MKL_LIB ?=


CC=icx
CFLAGS=   -static  -Wall -O0   
CFLAGS+=$(if $(GSL_DIR),-I$(GSL_DIR))
CFLAGS+=$(if $(MKL_DIR),-I$(MKL_DIR))
CFLAGS+=$(if $(GSL_LIB),-L$(GSL_LIB))
CFLAGS+=$(if $(MKL_LIB),-L$(MKL_LIB))

CLINK=  -qmkl   -lgsl -lgslcblas -lm

EXTFLAG = -g -fsanitize=address 

SRC_DIR=src
OBJ_DIR=obj

SRCS=$(wildcard $(SRC_DIR)/*.c)
OBJS=$(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
$(shell mkdir -p $(OBJ_DIR))


targets=MyTm


$(targets) : $(OBJS)
	@$(CC) $(CFLAGS)   -o ./$@  $^ $(CLINK)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	@$(CC) $(CFLAGS)  -c $< -o $@

clean :
	rm -rf $(OBJ_DIR) $targets

