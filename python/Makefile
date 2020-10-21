include=$(patsubst %, -I$(root_dir)/src/%, $(src_dir))
include2=$(patsubst %, /I$(root_dir)/src/%, $(src_dir))
vpath %.exe $(root_dir)
vpath %.o $(obj_dir)
All:Euler1D_wrap.o clean
Euler1D_wrap.o: Euler1D_wrap.cpp
	cl.exe  $(SOURCES) Euler1D_wrap.cpp  /ID:\\Anaconda32\\envs\\spam\\include /I"C:\\Program Files (x86)\\Windows Kits\\10\\Include\\10.0.18362.0\\ucrt" \
                /I"C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.27.29110\\include"\
				/I"C:\\Program Files (x86)\\Windows Kits\\10\\Include\\10.0.18362.0\\shared"  $(include2) \
				/link "/LIBPATH:C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.27.29110\\lib\\x64"\
                "/LIBPATH:D:\\Anaconda32\\envs\\spam\\libs"  "/LIBPATH:C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Tools\\MSVC\\14.27.29110\\lib\\x64"\
                "/LIBPATH:C:\\Program Files (x86)\\Windows Kits\\10\\Lib\\10.0.18362.0\\um\\x64"  "/LIBPATH:C:\\Program Files (x86)\\Windows Kits\\10\\Lib\\10.0.18362.0\\ucrt\\x64" \
				/LiBPATH:$(obj_dir) /DLL /OUT:$(root_dir)\_Euler1D.pyd
Euler1D_wrap.cpp: 
	swig  -python -c++ -o Euler1D_wrap.cpp $(include)  Euler1D.i
.PHONEY: clean
clean: 
	del *.obj
	del *_wrap.cpp
	del $(root_dir)\*.lib
	del $(root_dir)\*.exp
	move Euler1D.py $(root_dir)