objs=$(wildcard *.o)
vpath %.o $(obj_dir)
vpath %.exe $(bin_dir)
main.exe: 
	$(CPP_C) -o $(root_dir)/$@ $(objs) $(CPP_FLAGS)
