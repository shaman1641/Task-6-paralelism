all:cpu gpu
	
gpu: main_gpu.cpp
	pgc++ -acc=gpu -Minfo=all -lboost_program_options -o gpu main_gpu.cpp

cpu:serial multi

serial: main_cpu.cpp
	pgc++ -acc=host -Minfo=all -lboost_program_options -o serial main_cpu.cpp

multi: main_cpu.cpp
	pgc++ -acc=multicore -Minfo=all -lboost_program_options -o multi main_cpu.cpp

clean:all
	rm serial multi gpu