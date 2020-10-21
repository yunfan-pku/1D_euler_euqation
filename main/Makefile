include=$(patsubst %, -I$(root_dir)/src/%, $(src_dir))
vpath %.exe $(root_dir)
vpath %.o $(obj_dir)
main.o: main.cpp
	$(CPP_C) main.cpp -c $(CPP_FLAGS) -o $(obj_dir)\$@ $(include) 