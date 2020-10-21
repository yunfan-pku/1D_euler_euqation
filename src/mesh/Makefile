object = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
include=$(patsubst %, -I$(root_dir)/src/%, $(src_dir))
vpath %.o $(obj_dir)
vpath %.exe $(root_dir)
$(object): %.o:%.cpp %.hpp
	$(CPP_C) -c $< -o $(obj_dir)/$@ $(include) $(CPP_FLAGS)